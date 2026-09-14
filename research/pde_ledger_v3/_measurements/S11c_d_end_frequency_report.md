# S11c-d end-frequency checkpoint

## Current repair refresh

Both end-frequency packets have been recomputed from the full native producer
using the repaired c2 trace binding. The actual operator/elimination fingerprints,
isolated-root/lift records and source/input pins join at both ends. Each packet
contains 18 candidates, 22 full basis directions, 18 full-rank frequency pairings
and 45 exact-zero reconstruction/certificate scalars. Each end still has four
exactly real and fourteen nonreal normal candidates at the supplied binding.

| Refresh result | LEFT | RIGHT |
| --- | ---: | ---: |
| Wall seconds / peak RSS KiB | 30.73 / 201012 | 49.84 / 201044 |
| Largest algebraic coefficient-frame residual norm | 1.968e-13 | 2.103e-13 |
| Frequency-difference error, relative step 1e-4 | 2.564e-7 | 2.544e-7 |
| Frequency-difference error, relative step 5e-5 | 6.397e-8 | 6.356e-8 |

Both publication validators reconstruct 508 objects, 1,018 tags and 6,573
metadata paths per end. The finite-step residuals above remain nonzero.
The independent exact normal-reality certificate and complete subspace/rank
checks remain in each packet. These are pointwise algebraic frequency data,
not physical-current normalization or section 3b frequency-pole/Riesz data.

- [LEFT refresh](../scripts/out/S11c_d_end_frequency_left_trace_repair.out),
  793,815 bytes; [inventory](S11c_d_end_frequency_left_trace_repair_checkpoint.json).
- [RIGHT refresh](../scripts/out/S11c_d_end_frequency_right_trace_repair.out),
  793,349 bytes; [inventory](S11c_d_end_frequency_right_trace_repair_checkpoint.json).

The [c2 repair](S11c_c2_trace_repair_report.md) has cleared the source discrepancy
at both ends; LEFT also preserves all 85 baseline source objects. The next
implementation is the [both-end two-frequency physical current](S11c_d_both_end_current_plan.md),
followed by variable-profile matching. Earlier packets and the report below
retain their original input/provenance and are historical.

## Earlier end-frequency checkpoint

Both LEFT and RIGHT LAB_HELD/RHO4_CONSTANT frequency calculations are complete
at the supplied profile/parameter binding. The separate right acoustic-current
source calculation was unfinished when these frequency packets were built. It
has since [completed on this box](S11c_d_end_current_source_report.md), exposing
ten retained first-contrast source discrepancies. These frequency calculations
do not depend on resolving that physical-current source join.

`EndModeFrequencyData` extends the existing native engine (class line 5132;
construction line 5159; emission line 5297). All earlier native source nodes,
including the import wiring, Fourier reduction and current constructors, are
unchanged from the preceding source checkpoint. The focused runner consumes
the repaired full producer's original rational end pencils, root disks and
both momentum lifts. Source-method, input, profile-limit, bound-pencil and
elimination fingerprint joins hold at both ends. It does not rerun isolation
or replace the full eigenspace by a single vector.

Frequency stays live through the material binding and differentiation. The
computed radical transport supplies the total frequency derivative of the
nonlinear pencil. At every native root/lift the construction computes the
complete left/right SVD subspaces, full frequency pairing, normalized left
basis, and frequency-weighted algebraic field projector. Every degenerate
space also undergoes nonunitary changes of both full bases, with pairing,
normalization and projector covariance residuals retained.

| Measured result | LEFT | RIGHT |
| --- | ---: | ---: |
| Isolated radical disks / normal lifts | 9 / 18 | 9 / 18 |
| Scalar / two-dimensional subspaces | 14 / 4 | 14 / 4 |
| Complete basis directions | 22 | 22 |
| Full-rank frequency pairings | 18 | 18 |
| Exactly real / nonreal normal candidates | 4 / 14 | 4 / 14 |
| Exact scalar reconstruction residuals, all zero | 45 | 45 |
| Wall seconds / peak RSS KiB | 31.67 / 200692 | 66.17 / 200652 |

There are **3,920 numerical algebraic residual scalars** across both ends:
full kernels, nonlinear normalization, projector rank/range/idempotency and
basis covariance. The largest coefficient-frame residual norm is
**2.535e-13**. A further **1,800 derivative-check scalars** compare against
centered frequency differences transported on the same local radical branch.
The maximum errors at relative steps 1e-4 and 5e-5 are respectively
2.564e-7 and 6.397e-8 (LEFT), and 2.547e-7 and 6.376e-8 (RIGHT).
These finite-step residuals remain nonzero in the transcript.

The existing exact reality algorithm reconstructs the actual square-free
factors, joins axis intervals to every isolated disk, excludes denominator
intersections at this binding, and reports the normal-square sign. Independent
certificate validation recounts those roots and checks interval containment
and sign bounds. The original sheet records remain attached per candidate.
These are checks at the supplied binding, not a constant-rank theorem over
parameters or global sheet/exceptional-locus coverage.

The algebraic field projector is **not** a computed contour Riesz projector
or section 3b frequency-pole result. No physical-current normalization or
current-based incoming/outgoing assignment is made here. Real normal momentum
alone is not promoted to an open flux channel. The complete scattering,
profile-frequency bound-pole search, bookkeeping, controls and export remain
unfinished under the retained contract and supplied upstream debts.

## Published artifacts and validation

Each run exited 0 with empty stderr. Each transcript contains **508 computed
objects, 1,018 tags, 6,573 dimension/grade paths and 508 fresh lowerCamel keys**.
Validation reconstructs every literal or SHA/PIT payload, the complete
dimension/multigrade/lambda metadata and source-line index. It separately
rebuilds the frequency pairing/projector and checks full subspace ranks.

- [LEFT transcript](../scripts/out/S11c_d_end_frequency_left.out), 793,817 bytes;
  [source/run/publication inventory](S11c_d_end_frequency_left_checkpoint.json).
- [RIGHT transcript](../scripts/out/S11c_d_end_frequency_right.out), 798,813 bytes;
  [source/run/publication inventory](S11c_d_end_frequency_right_checkpoint.json).

Publication used fresh scratch outputs and destination-local atomic replacement;
earlier annexed outputs were preserved. These new outputs are uncommitted and
must use DataLad/git-annex at a future requested checkpoint. The one development
case has both end packets; no four-case production regeneration or export was
performed. No review/comparator/Wolfram/downstream run or authority/S10/Lean
edit occurred. The retained solver/export contract is byte-identical.

## Runnable checkpoint and next work

From the ledger directory, with the pinned producer caches available:

```bash
s11cdFrequencyRun=$(mktemp -d /tmp/s11cd-end-frequency-rerun.XXXXXX)
python -u _measurements/S11c_d_end_frequency_check.py \
  --manifest /tmp/s11c-mechanical-repair-20260912/d_full/manifest.json \
  --input _measurements/S11c_d_channel_preflight_input.json \
  --end RIGHT --run-directory "$s11cdFrequencyRun" \
  > "$s11cdFrequencyRun/full.out" 2> "$s11cdFrequencyRun/stderr.txt"
python _measurements/S11c_d_end_frequency_validate.py \
  --run-directory "$s11cdFrequencyRun"
```

Use `--end LEFT` for the other end. Publication refuses an existing destination.
Next prepare the variable-profile operator and unnormalized matching assembly
from the actual reduced rows and full vertex. Retain both-end and evanescent
mode data, explicit nonlocal integration/tail domains and refinement checks.
Keep physical flux outputs conditional on resolving the current-source discrepancy;
no computed current or pole result may be inferred from this preparation.
