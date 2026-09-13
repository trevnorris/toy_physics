# S11c-d two-frequency current checkpoint

The approved acoustic/face-lift split is implemented in the native
`ClosedCurrentPairing` builder. The focused checkpoint has **802 literal zero
scalar residuals in 34 families**, with no unresolved dimension constraints.
The former bulk cancellation stall does not require another upstream repair
on the evidence obtained here. Repair checkpoint `c643112a` remains committed;
this continuation is uncommitted.

| Computed residual group | Scalars |
| --- | ---: |
| Local acoustic polarization and both wave-curve divisions | 7 |
| Full face rows, linear lifts, matrix reconstructions and port joins | 534 |
| Finite-depth compositions, denominator numerator and boundary identities | 152 |
| Slab balance and same-frequency current join | 50 |
| Regular-sheet derivative checks and equal-depth integral derivative | 59 |

The local acoustic proof is reduced on the two independently computed wave
equations before closing the face amplitudes. The equal-depth-momentum branch
uses a separate wave eliminant. Each face constitutive check covers the whole
five-field linear row; the bulk and port matrices are separately reconstructed
from those rows. The finite-depth balance combines these certificates with
the slab balance, retaining interface power and the upper-boundary current.
Its unequal-depth composition is checked by clearing the actual inverse
carriers and expanding the numerator. The original symbolic denominator is
emitted; cancellation does not extend its domain.

The local acoustic identities and structural lift/port reconstructions retain
symbolic material parameters. Slab balance, face constitutive row reductions
and the final composed balances use the explicit LAB_HELD/RHO4_CONSTANT
reference material binding, leaving both frequencies and all four normal/bulk
momenta live. This is an exact identity check on that slice, not arbitrary
material-parameter coverage. The input has inactive V/X drivers and equal
memory times; the constructed source operands retain the independent kernels.

Normal-momentum and frequency derivatives are computed with the acoustic
radical chain rule. Their source operands retain all four product-rule terms;
the balance operands retain energy, current, interface exchange and the top
boundary. The local differentiated wave division and source-product
reconstruction are checked. The equal-depth first derivative of the integral
is computed separately. These records do **not** complete the mode-contracted
current/frequency-pairing relation or flux normalization. The radical derivative
requires its emitted denominator to be nonzero. No threshold/defective-mode
coverage or global sheet atlas is inferred. Infinite-depth use still requires
the previously derived decay domain; the equal-depth finite integral is not
an infinite-depth prescription.

The completed calculation took 221.17 seconds and 126,760 KiB peak RSS. Output
validation then found three missing unit tags on zero derivative operands.
After those tags were repaired, a source-checked reuse of the calculation
produced the guarded publication in 65.26 seconds and 123,472 KiB peak RSS,
with exit code zero and empty stderr. Cache reuse checks constructor and proof
method ASTs, input/source pins and artifact hashes. The publication has 234
tags (117 object/metadata pairs), 232 indexed source assignments, complete
metadata and a checked fingerprint/literal-residual inventory.

The atomically replaced focused transcript is
`scripts/out/S11c_d_modal_current_check.out`: 1,386,029 bytes, SHA-256
`9a2cee4e624b0be382931588ea51bdaa1a155744755922e43c7b0b5210a1fed5`.
The [checkpoint inventory](S11c_d_modal_current_checkpoint.json) records the
sources, proofs and publication. Heavy objects and report hashes use DAG-based
carrier/PIT fingerprints, avoiding full symbolic expression serialization.
The transcript is to be annexed at the next requested commit.

Next is full-subspace mode contraction and the current/frequency-pairing
relation, followed by conditional flux normalization and source controls; see
[the continuation plan](S11c_d_modal_current_plan.md). The committed four-case
production transcript remains the repair checkpoint. All ten broad engine
TODOs remain; no incomplete export or new four-case run was produced. Supplied
upstream debts remain explicit. Constant-end momentum poles remain distinct
from §3b's profile-dependent frequency poles. No S10/Lean or authority edits,
review legs, comparator, Wolfram, downstream steps, commits or pushes occurred.

## Runnable guarded checkpoint

Run from the ledger directory into a fresh destination. Never redirect
through an annex link. This command recomputes the new pairing and proofs;
the stable repaired current is source-checked before reuse.

```bash
mkdir -p /tmp/s11cd-current-rerun
python -u _measurements/S11c_d_modal_current_check.py \
  --manifest /tmp/s11c-mechanical-repair-20260912/d_reference/manifest.json \
  --input _measurements/S11c_d_channel_preflight_input.json \
  --current-manifest /tmp/s11c-mechanical-repair-20260912/d_current/manifest.json \
  --cache-result /tmp/s11cd-current-rerun/objects.pickle \
  --checks-json /tmp/s11cd-current-rerun/checks.json \
  --split-balances --derivatives --emit --require-zero-residuals \
  > /tmp/s11cd-current-rerun/full.out
```

The preserved publication is `/tmp/s11cd-split-current-20260912/validated`.
The full calculation is preserved in the sibling `published` directory;
its initial transcript failed metadata publication validation, not the
calculated residual checks. Optional `--construction-cache` and
`--verification-cache` arguments reuse only compatible, hashed computations.
Omit `--current-manifest` to reconstruct the earlier current stage as well.
