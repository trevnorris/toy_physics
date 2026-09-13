# S11c-d reference mode-subspace current checkpoint

Checkpoint `5b007ddb` saved the preceding two-frequency current proof, using
Git for source/reports and DataLad/git-annex for its transcript. This subsequent
reference-subspace continuation is committed as `7ba922a9`. No upstream repair or new physical premise is
indicated by its computed residuals.

The native `ModalCurrentSubspaces` builder reconstructs **all 18 root/lift
candidates** in the repaired LAB_HELD/RHO4_CONSTANT reference packet: 14 scalar
and four two-dimensional nullspaces. Every right and left basis has the full
producer nullity, with basis, kernel and coordinate-projector residuals.
Every frequency-pairing matrix has full computed rank. This continues the
isolated-root checks; it does not replace them with a generic witness.

`prepare` differentiates the actual two-frequency current, source and closed
pencil operands along the acoustic wave curve, including the radical and
depth-integral derivatives. `exact_low_degree_lifts` factors the bound pencil
and joins its degree-two factor to the existing disjoint disks. Four candidates
obtain exact normal/bulk momenta and zero factor/wave residuals. The other
14 retain their isolated numerical roots. No new input was chosen to create
channels. The radical transport is restricted to its nonzero denominator.

`construct` forms the complete matrices `Nω = L† Pω R` and `Nk = L† Pk R`.
It independently contracts the physical S11b slab-plus-bulk current and energy
on the right-field basis. The actual plus-row power map gives the covectors
`C = R† B+`. Their computed decomposition into `G L† + D` retains the full
projection defect `D`. The differentiated balance then reconstructs the
physical current/energy from the pencil pairings and the remaining source,
interface, top-boundary and phase-derivative terms. The split residual is an
algebraic reconstruction check; the balance is supported by the preceding
802-residual current proof. A row-dual left vector is not silently interpreted
as a physical test field.

| Numerical residual family | Maximum matrix norm |
| --- | ---: |
| Right / left kernel | 5.50e-15 / 4.56e-15 |
| Finite-depth balance | 6.03e-15 |
| Normal / frequency balance derivative | 8.01e-13 / 2.15e-14 |
| Reconstructed physical current / energy | 8.21e-13 / 2.20e-14 |
| Frequency normalization | 5.56e-15 |
| Signed-current / flux-frequency normalization | 3.24e-16 / 2.69e-16 |

These rounded-up norms are double-precision implementation diagnostics in
the supplied reference-unit coordinates, not the withheld physical conversion
criterion. SVD nullity uses `1e-8 max(1,s_max)`; frequency/current ranks use
`1e-9 max(1,norm_2)`. The inventory retains all 27 residual families and their
2,794 scalar entries. The 250 quadratic-extraction scalar residuals and the
exact spectrum reconstruction residuals are zero.

Finite-depth checks use the input `W_0 = 1` as a diagnostic cutoff. It is not
an additional physical boundary. Infinite-depth integrals are evaluated only
for the eight candidates with disk-certified positive bulk decay. Full current
matrices and their anti-Hermitian residuals are retained before diagonalization.
Two physical-sheet, exactly real-normal, two-dimensional spaces (records 16
and 17) admit signed physical right-current normalization. Their unnormalized
current eigenvalues are approximately `+0.58896094947` and `-0.58896094947`,
each twice, in the declared unit frame. The resulting four basis directions
have separately emitted field/flux maps, left/right kernel residuals and full
frequency-normalization residuals. Growing, nonphysical-sheet and complex-normal
candidates remain matching data. No end-oriented incoming/outgoing assignment
or sector classification is inferred here.

The epsilon-squared physical forms and their coefficient-based normalization
maps are distinguished in the metadata. Eigenvectors restore their field or
dual-row units; numeric SVD/projector diagnostics use the declared coordinate
frame and do not define an invariant physical Euclidean metric. Recovered
polynomial factors and their zero residuals restore frequency-squared units.
All values are evaluated at the reference `eta = sigma_W = 0`; independent
material-response kernels remain in the preceding symbolic construction.

The final run took 128.69 seconds and 145,324 KiB peak RSS, with exit code zero
and empty stderr. A unit-label correction rerun reproduced every operand,
form and residual tensor exactly. The checked transcript has 2,660 tags,
1,330 object/metadata pairs, 2,658 source assignments and 1,298 tensors checked
against the retained payload, including SHA/PIT or literal comparisons and
10,528 tensor metadata paths. The
[checkpoint inventory](S11c_d_modal_subspace_checkpoint.json) records source,
input, cache and transcript hashes. Publication is atomic under
`scripts/out/S11c_d_modal_subspace_check.out`; it was annexed in `7ba922a9`.
The original current and four-case production transcripts remain
the committed checkpoints.

Next: source-level mass-transfer/bulk-load/boundary-work controls, the canonical
adjoint-row/physical-field current export join, then both ends and the other
cases; see the [continuation plan](S11c_d_modal_current_plan.md). This does not
complete the mode/flux TODO, global sheet or exceptional-domain coverage,
variable-profile matching, scattering, profile-frequency poles or export.
Constant-end normal-momentum roots are distinct from section 3b's conditional
profile-frequency bound poles. Supplied section 1 premises and upstream debts
remain supplied. No S10/Lean or authority edit, review leg, comparator, Wolfram,
downstream run or push occurred.

The later source controls exposed a higher-degree reality-certificate gap.
The approved repair is now complete; see the
[source-control report](S11c_d_current_source_controls_report.md) and
[local repair plan](S11c_d_current_reality_repair_plan.md). A new reference rerun
certifies all 18 normal-reality statuses while preserving all 1,298
earlier tensors and the four normalized basis directions. Its separate
[inventory](S11c_d_modal_reality_reference_checkpoint.json) and transcript
preserve this original committed checkpoint.

## Runnable focused checkpoint

Run from the ledger directory into a fresh scratch directory, preserving the
referenced source-pinned caches. Never redirect through an annex link.

```bash
mkdir -p /tmp/s11cd-subspace-rerun
python -u _measurements/S11c_d_modal_subspace_check.py \
  --manifest /tmp/s11c-mechanical-repair-20260912/d_reference/manifest.json \
  --input _measurements/S11c_d_channel_preflight_input.json \
  --current-manifest /tmp/s11c-mechanical-repair-20260912/d_current/manifest.json \
  --current-checkpoint _measurements/S11c_d_modal_current_checkpoint.json \
  --spectrum-manifest /tmp/s11c-mechanical-repair-20260912/d_full/manifest.json \
  --run-directory /tmp/s11cd-subspace-rerun --emit \
  > /tmp/s11cd-subspace-rerun/full.out
```

The preserved run is `/tmp/s11cd-modal-subspaces-20260912/validated`. Its manifest
records the separate post-emission publication checks; the command above
recomputes and emits the objects and residuals without publishing over a
tracked output.
