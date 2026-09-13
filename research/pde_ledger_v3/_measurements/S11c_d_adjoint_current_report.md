# S11c-d adjoint-row / physical-field current map

Checkpoint `c7f2d879` commits the exact normal-reality repair and completed
source controls. Its three `.out` files are DataLad/git-annex objects; Git
stores their pointers. The continuation below is new, uncommitted builder
work on the same LAB_HELD/RHO4_CONSTANT reference input.

`AdjointCurrentMap` now constructs the explicit field representation of the
frequency-normalized row dual on **all 18 isolated root/lift candidates**.
It retains all **22 basis directions**: 14 scalar spaces and four complete
two-dimensional spaces. No input, upstream authority or physical premise
was changed. This completes the focused regular-reference mapping step;
both-end current construction and profile matching remain next.

## Computed map and current reconstruction

Let P be the physical pencil, B the plus-row power map derived from the
source-work ansatz, R the full right basis, and L the existing
frequency-normalized row-dual basis. The construction solves
`B† A = L` and forms `Q = B P`. Thus A is an adjoint field representation
for the power-weighted row representation Q, with the emitted residuals
`Q† A` and `Q R`. It is not identified with an untransformed adjoint of P.
Row maps and derivative pairings use the established epsilon-squared-removed
coefficient convention; physical current/work forms restore epsilon squared.

Both derivatives are the existing **right-leg partial/radical derivatives**,
holding the independent left leg fixed. Direct differentiation of Q agrees
with `Bν P + B Pν`; the `Bν P` term is retained in the modal pairing.
The 50 symbolic product-rule residual scalars are exactly zero. The finite-depth
mixed current and energy are also reconstructed from the full differentiated
source/balance identity, including interface, upper-boundary, phase, depth
integral and other source-work terms. This uses the previously verified
[two-frequency current proof](S11c_d_modal_current_checkpoint.json).

The existing covector decomposition `R† B = G L† + D` supplies the computed
bridge G and projection defect D. Solving `B† F = D†` gives
`R = A G† + F`. The physical forms are reconstructed as
`R† J R = G (A† J R) + F† J R`, with the analogous energy identity.
Every term is retained. Eight disk-certified decaying candidates also receive
the infinite-depth current forms. The two originally eligible physical-sheet
spaces retain their signed-current normalization: four directions.

The distinction between mixed and physical current is numerically visible:
on reference records 16/17, the mixed current is approximately
`±1.17792189894 i I₂`, while the bridge is approximately `−0.5 i I₂` and
the physical current is `±0.58896094947 I₂`, in their respective restored
reference-unit coefficients. The complex phase is part of the computed map.
No phase or projection defect is discarded to identify these different forms.

A nonunitary basis change is applied to every complete subspace. Recomputed
maps verify the similarity transformation of the mixed current, congruence
transformation of physical current, bridge/defect covariance and frequency
normalization. This is a coordinate regression of the computed construction.

## Validation and publication

| Quantity | Measured result |
|---|---:|
| Defined field maps / candidates | 18 / 18 |
| Power-map rank | 5 on each candidate |
| Complete adjoint-field basis directions | 22 |
| Exact symbolic product-rule zeros | 50 |
| Literal numerical residual scalars / families | 3,056 / 29 |
| Largest numerical residual norm | 3.679e-13 |
| Physical-field reconstruction maximum | 3.576e-16 |
| Signed-current reconstruction maximum | 1.666e-16 |
| Verified tensor payloads / metadata paths | 1,410 / 11,505 |
| Output tags / metadata pairs | 2,876 / 1,438 |

Power-map rank uses `1e-10 max(1, largest singular value)` in the declared
reference-unit coefficient frame. The publication check verifies adjoint-field
rank against each full source nullity at tolerance `1e-9`; the numerical
reconstruction diagnostic bound is `1e-8`. These are pointwise numerical
checks, not a constant-rank theorem over parameter space. The prior exact
normal-reality, sheet/decay and normalization domains remain intact.

The run exited 0 with empty stderr in **157.53 seconds**, using **125,480 KiB**
peak RSS. The publication validator checked source/input/cache pins, every
root/lift join, symbolic carrier fingerprints, numerical SHA and three
position-weighted PIT samples per heavy tensor, literal residuals, restored
units, metadata paths and the final emission index. There are no unresolved
emitted dimension constraints. Full cached objects were saved before emission.
Only one heavy CAS process ran at a time.

The transcript was published through a destination-local temporary file and
atomic replacement at
[`scripts/out/S11c_d_adjoint_current_check.out`](../scripts/out/S11c_d_adjoint_current_check.out):
**1,206,643 bytes**, SHA256
`1c248e5f16f184ff5cef987289a5ece0f0e6f4a4dc89522f31647a8be13a7abd`.
The [inventory](S11c_d_adjoint_current_checkpoint.json) records source hashes,
payload pins, all numerical residual maxima and per-mode map summaries.
The native change adds only `AdjointCurrentMap`; earlier native definitions
are source-checked against the committed modal producer. Committed transcripts
are preserved. The new transcript awaits annex storage at the next requested
checkpoint.

## Remaining work

Extend the verified full-subspace current and field maps to both asymptotic
ends of the development case, preserving all root, sheet and rank outcomes.
Then construct variable-profile matching from the retained operator and
re-expand the continuum response at the declared orders. Full-case integration,
global/exceptional-domain coverage, profile-dependent frequency bound poles,
controls, bookkeeping and the complete export remain open. Constant-end
normal-momentum poles do not answer the section 3b bound-pole question.
All ten broad native TODOs remain; this is not engine completion.

No S10/Lean/authority edits, review legs, comparator, Wolfram, downstream run,
full four-case regeneration or new export belongs to this focused result.
No further commit or push was made after `c7f2d879`.

## Reproduce this checkpoint

From the ledger directory, the following uses the pinned upstream scratch
caches and a fresh run directory. It validates a scratch transcript; publication
is a separate guarded operation that refuses an existing destination.

```bash
s11cdRun=$(mktemp -d /tmp/s11cd-adjoint-rerun.XXXXXX)
mkdir -p "$s11cdRun/source"
cp --parents scripts/S11c_d_mixing_scattering_sympy_audit.py \
  _measurements/S11c_d_adjoint_current_check.py "$s11cdRun/source/"
python -u _measurements/S11c_d_adjoint_current_check.py \
  --manifest /tmp/s11c-mechanical-repair-20260912/d_reference/manifest.json \
  --input _measurements/S11c_d_channel_preflight_input.json \
  --current-manifest /tmp/s11c-mechanical-repair-20260912/d_current/manifest.json \
  --current-checkpoint _measurements/S11c_d_modal_current_checkpoint.json \
  --modal-checkpoint _measurements/S11c_d_modal_reality_reference_checkpoint.json \
  --run-directory "$s11cdRun" --emit \
  > "$s11cdRun/full.out" 2> "$s11cdRun/stderr.txt"
python _measurements/S11c_d_adjoint_current_validate.py \
  --run-directory "$s11cdRun" \
  --modal-checkpoint _measurements/S11c_d_modal_reality_reference_checkpoint.json
```
