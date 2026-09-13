# S11c-d right-current runtime investigation — 2026-09-13

The user clarified that a run lasting several hours is acceptable if this
box can complete it. Reactivate the required right-end source construction;
the previous deferral and stop inventories remain historical evidence.

1. Isolate the actual mechanical closure increment from the pinned RIGHT
   LAB_HELD/RHO4_CONSTANT reduced pencil. Save its unsimplified operands and
   input/source hashes before testing simplification routes. Measure CPU,
   wall time and peak RSS; inspect the stage that is actually running.
2. Investigate algebraically equivalent evaluation orders first, keeping
   all material parameters, frequency, momentum and background grades live.
   Preserve every input denominator exclusion. Do not drop higher grades or
   interpret an interrupted cancellation as a residual result. Exact
   coefficient extraction must not silently invoke the original costly
   all-variable cancellation itself.
3. Add durable stage/entry checkpoints to the right-end run. Resume only when
   source/input/operand hashes match. Permit an hours-long serial run with
   periodic progress/resource records and stack diagnostics. Never launch
   two heavy CAS processes concurrently or redirect through annex links.
4. If an equivalent algebraic refactor is used, verify reconstruction from
   the original unsimplified numerator/denominator operands, retain excluded
   loci, and compare the completed LEFT construction with its prior checkpoint.
   Finish the RIGHT conservative, mass-rate and acoustic source objects;
   emit their actual raw and retained-grade residuals before any guard.
5. Validate and atomically publish the bounded right source transcript and
   record the runtime/resources. Update the deferred tracker and current
   plans according to the computed result. Continue the physical current
   extension only if the required source joins are established. A nonzero
   retained-grade source discrepancy requires a separate physics assessment.

No new physical premise, profile tuning, upstream/authority edit, review leg,
comparator, Wolfram run, downstream step or commit. The recent both-end frequency
packets and the retained solver/export contract remain preserved.

## Execution notes

The initial isolated mechanical increment completed in 26.79 seconds while
keeping all material symbols. The historical stop's source line 2175 actually
names the mass increment; its inventory's mechanical label was one line off.
A monitored full attempt also found dense cancellation in the face load.
The exact collector now separates field amplitudes and background powers,
retains background-dependent denominator factors explicitly, and applies the
material cancellation to each resulting coefficient. It performs no background
truncation; the original retained projections remain at their source sites.

The expanded full attempt saved 125 rational operations, then Python segfaulted
in a recursive dense GCD for the final grade of the mechanical-sum diagnostic.
The observed process peak was about 201180 KiB, with an 8 MiB stack limit.
The 256 MiB stack retry completed all 127 rational operations and saved the
full acoustic construction. It then segfaulted at a different site: eager
expression expansion in the first equivalence check. The stack hypothesis is
not established. A direct sparse-ring conversion completed that same identity
in 30.37 seconds with a zero cross-product residual. The final run preserves
the completed native diagnostic sum, seeds only hash-verified construction
packets under identical native/loader/input pins, and recomputes every
equivalence check with direct sparse-ring conversion. Frozen sources, failure
traces, raw operands and per-operation hashes remain under
`/tmp/s11cd-right-source-resume-20260913`. No material binding was needed.

The first sparse verification run saved 96 exact-zero identities, then spent
several minutes expanding the last mechanical-increment reconstruction. A
bounded structure inspection found numerator/denominator total degrees
154/146 for the uncombined source fraction versus 21/13 after `together`;
that exact fraction combination took 0.66 seconds. The run was interrupted
with exit 130, retaining its saved identities. Verification now combines
fractions before direct sparse-ring conversion. It emits the original
before/after expressions and their uncombined denominator exclusions as well
as the combined numerator/denominator operands and cross-product residual.
The construction packets remain unchanged; every proof is recomputed under
the new instrument pin. No material substitution or reduced-order shortcut
is introduced by this verification change.

All 127 combined-fraction identities computed exact-zero residuals. Emission
then exposed a metadata boundary: an untruncated rational background expression
cannot be assigned finite polynomial support. Its metadata now records the
computed exact numerator/denominator supports and the full fraction's restored
unit, with no invented finite-order truncation. Zero residual units are derived
from the physical source rows and the actual denominator factors. The final
serialization run may seed proof packets only if the exact identity-function
AST is unchanged, in addition to the existing construction/input/operand pins.


## Completed checkpoint

Both source packets and their complete arithmetic proofs are emitted and
validated. RIGHT: 127 zero arithmetic identities; ten nonzero source-join
entries at `(0,1,0)`. LEFT: 82 zero arithmetic identities and zero source joins;
all 85 source objects agree with the baseline. Parameter alignment is compared
by its full symbolic key set, with duplicate-key rejection, rather than by
incidental atom-set iteration order. Original transcripts remain preserved.
The source report and runtime inventory record timings, source pins, failures,
publication checks and the physical-source boundary that remains open.

A fresh serial source run from the ledger directory can be reproduced with:

```bash
s11cdRuntimeRun=$(mktemp -d /tmp/s11cd-source-runtime.XXXXXX)
ulimit -s 262144
python -u _measurements/S11c_d_end_current_resumable.py \
  --manifest /tmp/s11c-mechanical-repair-20260912/d_full/manifest.json \
  --input _measurements/S11c_d_channel_preflight_input.json \
  --end RIGHT --run-directory "$s11cdRuntimeRun" \
  > "$s11cdRuntimeRun/full.out" 2> "$s11cdRuntimeRun/stderr.txt"
python _measurements/S11c_d_end_current_source_validate.py \
  --run-directory "$s11cdRuntimeRun"
```

RIGHT's computed source residual currently triggers exit 1 only after the
complete transcript and checks are written. `--resume` requires an unchanged
instrument/source/input signature; source-checked `--seed-run` records its
reuse provenance and checks exact proof-function identity before proof reuse.
No new physical source formula was adopted. The next task is to locate the
retained source discrepancy and establish the proper repair scope.
