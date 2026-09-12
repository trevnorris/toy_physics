# S11c-d current preflight after the mechanical-load repair

The original five-coefficient mechanical-load mismatch is resolved on the
freshly reduced LAB_HELD/RHO4_CONSTANT reference. All five mechanical
coefficients and all five mass coefficients agree. The existing d current
construction passes with the repaired b/c1/c2 exports; its formulas were
unchanged. The [repair record](S11c_mechanical_repair_report.md) documents the
upstream action-normalization correction and four-case checks.

[SlabEnergyBalance](../scripts/S11c_d_mixing_scattering_sympy_audit.py#L1870)
varies the inherited, tangentially reduced stored energy with all five fields
independent, then applies the instantaneous material virtual constraint. It
reconstructs the conservative variation and boundary work, retains the actual
density-rate defect, and computes its additional boundary current and chemical
derivative. All eleven variation, integration-by-parts, current, and polarization
residuals remain zero. The density-gradient correction is retained before
imposing zero material mass rate.

[ClosedAcousticEnergy](../scripts/S11c_d_mixing_scattering_sympy_audit.py#L2049)
derives acoustic energy from the supplied wave equation and pressure-work flux.
It computes finite-depth integrals and the positive-decay infinite-depth limit.
Equal normal momenta have a separate integral equal to the depth cutoff; no
finite infinite-depth flux is assigned on the nondecaying diagonal. Both face
closures use the energy-derived chemical derivative with all three relaxation
times independent. The five acoustic algebraic residuals, five mass joins,
five mechanical joins, and four face-equation residuals are zero.

The independent stiffness anchor is still `W_0^2` in both the source thickness
row and the energy-derived row, with zero difference and restored dimension
`[2,0,0]`. The [current inventory](S11c_mechanical_repair_d_current_checks.json)
records all 30 zero construction/comparison scalars plus this separate anchor.
It finds no metadata gaps, nonfinite objects, unresolved dimensions, or
unassigned zero-map dimensions. The preceding
[reference inventory](S11c_mechanical_repair_d_reference_checks.json) checks
39 reconstruction scalar digests against the computed literal-zero digest.

The fresh reference reduction took 297.26 seconds with 1,714,256 KiB peak RSS.
The current run took 180.97 seconds with 122,636 KiB peak RSS, empty stderr,
stable source pins, and exit zero with `--require-zero-residuals` enabled.
The [published current transcript](../scripts/out/S11c_d_nonlocal_current_reference_preflight.out)
contains 146 unique tags and 1,735,781 bytes (SHA-256
`d61e60cd06dd1f39f0528ff609a4a224d0c278e4e3e79318047947ab63f543ae`).
[Publication provenance](S11c_mechanical_repair_d_publications.json) records
both fresh transcripts and preservation of the previous annex payload.

The earlier disagreement and its original input pins remain recorded in
[S11c_d_nonlocal_current_runs.json](S11c_d_nonlocal_current_runs.json) and
[S11c_d_nonlocal_current_inventory.json](S11c_d_nonlocal_current_inventory.json).
Those are historical pre-repair measurements; the published transcript now
contains the repaired result.

This is a reference-case current check. The affected four-case spectrum, sheet,
and threshold evidence has been regenerated, inventoried, and published; see the
[full repair inventory](S11c_mechanical_repair_d_full_checks.json). Variable-profile current/flux
normalization, complete two-ended scattering, profile-dependent frequency poles,
and the remaining engine/export program are still open. The supplied physics
and carried upstream debts are unchanged; these checks do not establish global
spectral or exceptional-locus coverage.
