# S11c-d current source controls — completed focused checkpoint

The approved local reality repair and all three source controls are complete
for the supplied LAB_HELD/RHO4_CONSTANT reference input. This is uncommitted
builder work following `7ba922a9`; it is not completion of S11c-d.

| Source control | Root/lift candidates | Exactly real / nonreal normal candidates | Normalized subspaces / basis directions | Largest numerical residual |
|---|---:|---:|---:|---:|
| `massTransferOff` | 18 | 8 / 10 | 4 / 6 | 2.622e-13 |
| `bulkDensityChange` | 18 | 4 / 14 | 2 / 4 | 2.019e-13 |
| `densityGradientOff` | 14 | 4 / 10 | 2 / 4 | 4.341e-13 |

Each candidate uses its entire computed left/right nullspace. Both momentum
lifts of every isolated radical root are retained. Current normalization also
requires the separate sheet, decay, nonzero-normal, denominator, Hermiticity,
frequency/normal pairing-rank and current-rank conditions. Ranks and matrix
residuals retain their stated numerical tolerances; normal reality is certified
by exact rational arithmetic on the isolated roots.

## Computation and controls

`CurrentSourceControls.build` re-enters at the energy/reduced-pencil sources.
The A/V control sets `Lambda_A_0 = Lambda_V_0 = 0`. The bulk-density control
uses a fresh positive density carrier bound to twice the supplied density; it
is an arithmetic/load-routing control, not the section 5 profile-FORM ablation.
The gradient control sets `kappa_theta = kappa_theta_W = 0`, retaining
`kappa_W`. It changes the characteristic degree and candidate count; the
inventory follows the actual polynomial rather than requiring 18 candidates.
This is finite-root coverage at that binding; no continuation through infinity
or coverage of all parameter strata is claimed.

Every control retains both source-energy and reduced-pencil operands, plus
base/control/change objects for 15 current/energy/power families. Five matrix sensitivities
use a common declared acoustic-wave point, explicitly separate from operator
modes. The mass-rate boundary correction is independently contracted on each
full mode space. See the [completed inventory](S11c_d_current_source_controls_completed.json)
for those literal differences, source substitutions, ranks, corrections and
all diagnostic maxima.

There are **2,259 literal zero source-reconstruction
scalars**, 753 per control in 25 families. The control maximum numerical
reconstruction diagnostic is 4.341e-13 in the reference unit frame.
The source computations, checks and isolated spectra of the two completed
packets were reused only after input/payload hashes, source substitutions,
bound-pencil equality and unchanged non-modal source AST checks. Reuse of
already repaired modal packets additionally checks the full modal/certificate
class ASTs and the exact extraction operands from their pinned transcript.
All earlier modal tensors in these packets compare bit for bit equal; only
the newly eligible mass-transfer-off flux tensors were added relative to the
original interrupted checkpoint.

## Reality repair and reference regression

`NormalRealityCoverage.construct` derives the radical-axis condition from the
actual bound wave relation. Per square-free factor it computes real/imaginary
coefficient gcds on both axes, isolates every real gcd root with rational
intervals, joins intervals strictly inside unique existing root disks, and
certifies the sign of the computed normal-momentum square. Off-axis disks
supply exact axis-exclusion evidence. Empty axis-root sets, the origin,
threshold zeros, denominator gcds and unresolved disk joins remain explicit.
Polynomials use the declared reference-unit coefficient frame; radical
intervals, disk margins and normal-square bounds restore their physical
units. Every emitted object carries dimensions and epsilon/eta/sigma/lambda
metadata. Symbolic carriers use SHA/PIT fingerprints; residuals are literal.

The A/V-off records 8, 9 now have computed flux maps.
Their current eigenvalues are 0.629795541945, -0.629795541945.
The opposite-sheet records remain excluded despite their certified real normal
momenta. The original reference was rebuilt against the same native spectrum:
all 1,298 earlier tensors are unchanged and its four
normalized basis directions remain. Across reference and controls,
**68 candidates** have explicit real/nonreal certificates with no unresolved
normal-reality status at these supplied bindings; **79 exact certificate
residuals** are zero.

The validation helper independently recounts interval roots, reconstructs
factor gcds, checks exact disk containment and encloses the normal-square
extrema. Algorithm-only fixtures deliberately exercise branch zero, threshold
zero, denominator intersection and an uncertified disk. They are not additional
physical-mode samples and do not establish global exceptional-locus coverage.
See the [reference inventory](S11c_d_modal_reality_reference_checkpoint.json)
and [certificate/fixture validation](S11c_d_current_reality_validation.json).

## Publication and scope

The complete control run returned 0 in 447.21 seconds, with
179,924 KiB peak RSS. Its stdout contains 8,252 tags,
4,126 metadata pairs and a verified final emission index.
Validation checked 3,614 modal tensors and
39,882 metadata paths, with no unresolved emitted dimension
constraints. The reference rerun returned 0 in 127.59 seconds
at 145,408 KiB peak RSS. Only one heavy CAS process ran at a time.

- `scripts/out/S11c_d_current_source_controls.out`: 5,959,452 bytes;
  SHA256 `6f70558bbf0b336050bf8acfa51a8ea0ad1e3127e1651500dc547e6ad9f60cac`.
- `scripts/out/S11c_d_modal_reality_reference_check.out`: 1,212,298 bytes;
  SHA256 `d5e1ae71c816bbd77894a90ee5696b0ae1a198474103cd424dbdf0d191f890a7`.

Both were published through destination-local temporary files and atomic
replacement. Existing annexed outputs were not redirected through or modified.
The original two-control interrupted transcript and its
[stop inventory](S11c_d_current_source_controls_checkpoint.json) are retained
as the diagnosis. Scratch attempts retain the repaired cache-guard, fingerprint
selection and cache-provenance serialization errors; none is presented as a
completed run.
No export or full four-case production transcript was regenerated at this
focused boundary. No S10/Lean/authority changes, reviewer, comparator, Wolfram,
downstream run, commit or push occurred.

Next: finish the canonical adjoint-row/physical-field current join with the
computed power map and projection defect, then extend to both ends and the
one-case profile-matching calculation. Global sheet/exceptional coverage,
section 3b profile-frequency bound poles, full controls/bookkeeping and the
complete export remain open under the retained solver contract.

## Runnable source-control checkpoint

From the ledger directory, use a fresh scratch directory. This command rebuilds
the controls from the pinned baseline; completed packets may optionally be
reused with `--reuse-checkpoint` and their inventory. The helper checks source
and payload compatibility before reuse.

```bash
mkdir -p /tmp/s11cd-controls-rerun
python -u _measurements/S11c_d_current_source_controls.py \
  --manifest /tmp/s11c-mechanical-repair-20260912/d_reference/manifest.json \
  --input _measurements/S11c_d_channel_preflight_input.json \
  --current-manifest /tmp/s11c-mechanical-repair-20260912/d_current/manifest.json \
  --current-checkpoint _measurements/S11c_d_modal_current_checkpoint.json \
  --run-directory /tmp/s11cd-controls-rerun --control all --emit --modes \
  > /tmp/s11cd-controls-rerun/full.out
```
