# S11c-d RIGHT source trace — c2 pressure trace/reference-slot mismatch

## Subsequent authorized repair

The [c2 trace repair](S11c_c2_trace_repair_report.md) is implemented and its
full four-case export/transcript is validated. Fresh symbolic LEFT/RIGHT source
checks now have zero retained residuals and all exact arithmetic identities
are zero; LEFT preserves all 85 baseline source objects. See the
[current source report](S11c_d_end_current_source_report.md). The investigation
below records the earlier input and diagnosis; its stop was superseded by the
user's authorization to repair. Full downstream d regeneration, its nine output inventories, both
end-frequency refreshes and the reference-current source refresh are complete
and published. The [both-end physical-current extension](S11c_d_both_end_current_plan.md)
is next.

## Historical investigation

Checkpoint `f55e55b6` commits all earlier S11 work, including the six annexed
transcripts and removal of the completed runtime item from the waiting list.
The subsequent source trace locates the ten RIGHT discrepancies in the
pressure-trace binding between c1 and the inherited slab row in c2. No physical
formula, upstream export or native d constructor has been changed.

The supplied c1 pressure is evaluated on the curved face, as specified in
`directives/S11c_c1_SHARED_PHYSICS.md:197` and `:241`. Its independent layer
construction includes the shifted pressure evaluation at
`scripts/S11c_c1_bulk_closure_sympy_audit.py:673`. In contrast, c-a constructs
`delta_p_plus/minus` and their normal jets as an affine bulk field around the
flat reference face (`scripts/S11c_a_interface_geometry_sympy_audit.py:578`),
then composes that field onto the physical face (`:588`, `:666`). Those source
slots and their shape terms reach the b slab rows.

In `scripts/S11c_c2_selfenergy_fold_sympy_audit.py:502`, `build_face` computes
the c1 response pressure. At `:505–513`, it anchors an outgoing extension of
that pressure at the reference face and substitutes both its value and normal
jet directly into the reference slots. The inherited face composition then
adds a displacement correction to a value that already represents pressure
on the curved face. This is the representation mismatch isolated here.

## Exact development-case evidence

The diagnostic works from the completed LEFT/RIGHT symbolic source packets and
their pinned b/c1/c2 exports. All material parameters remain symbolic. It
applies the supplied profile limits and differentiates the same harmonic field
ansatz; this is a provenance trace of the source rows, not a replacement of d's
required reduction by an unreduced operator.

- All five coefficients of the imported chemical operator agree exactly with
  the energy-derived chemical row at each end: ten zero residuals.
- The imported branch density agrees with the density extracted from the
  no-transfer mass rate at each end: two zero residuals. The RIGHT coefficient
  is `rho_br*(1+eta_bg)`; LEFT is `rho_br`.
- For each face, the diagnostic extracts the normal-pressure-jet coefficient
  directly from the imported b mass and mechanical rows. It computes the
  normal derivative of the outgoing continuation of the saved reference
  pressure, then contracts both faces and differentiates by all five amplitudes.
- The resulting inherited jet row equals the complete saved source discrepancy
  in both rows at both ends: twenty zero accounting residuals. RIGHT retains
  all ten original nonzero entries at `(epsilon,eta,sigma)=(0,1,0)`, lambda one;
  LEFT has zero jet contributions and zero source discrepancies.

Thus **32 exact scalar identities are zero**. In this development case the
additional pressure evaluation accounts for every observed discrepancy. This
does not establish the general mixed-grade repair or global parameter/sheet
coverage. No discrepancy has been subtracted from a native result to make a
current pass; both original rows remain unchanged and available.

The [instrument](S11c_d_right_source_trace_check.py) constructs the chemical
join at lines 104–109 and the pressure-jet accounting at lines 115–147. It
emits both operands, original discrepancies and literal residuals, with fresh
lowerCamel keys, restored units and multigrade/lambda metadata. Reciprocal
factors remain explicit domain operands; heavyweight values use native
carrier fingerprints and SHA. The [inventory](S11c_d_right_source_trace_checkpoint.json) records the
source/payload pins, all 50 fresh keys and validation census. The run exited zero
with empty stderr in **44.87 seconds**, at **1,620,708 KiB** peak RSS. Validation
reconstructed all **102 tags**, payload fingerprints, grade/unit metadata and
emission-line assignments, and checked the source pins again. The
[diagnostic transcript](../scripts/out/S11c_d_right_source_trace.out) has
**381,485 bytes**, SHA256
`81b8d1bd3776a12aa283aa8b6081cda0710096c9f697fd0deb15066f9845b04b`.
It was published atomically; prior transcripts retain their original hashes.

## Repair boundary

The next repair belongs at the **c2 source representation interface**, before
resuming d current normalization. Construct a map from c1's curved-face
pressure to the reference pressure/normal-jet slots expected by the inherited
trace ansatz, deriving it from that actual ansatz. Compose the map back through
the original trace to test it. Do not delete physical shape terms or patch the
d residual with its measured value.

First do the full retained `(eta,sigma)` rectangle, including mixed grade, in
one case; a constant-end correction alone is insufficient. Keep both momentum
legs, normal continuation, denominator domains, orientations and source order.
Then check both anchoring/density branches and regenerate the affected c2
exports/transcripts and d reduced spectra/current consumers. The current
frequency packets remain algebraic results for the old supplied pencil and
must be recomputed if that pencil changes. Present evidence does not require
changing the b energy, its mass density or the c1 physical face-pressure
contract; keep them fixed unless further source checks show otherwise.

Work stops at this upstream repair boundary in accordance with the user's
request to report a change in scope. The new diagnostic and plan are subsequent
to `f55e55b6`. No upstream repair, full production regeneration, new export,
review/comparator/Wolfram/downstream run or push has been performed.

## Reproduce

From the ledger directory, use the pinned end packets and a fresh directory:

```bash
s11cdTraceRun=$(mktemp -d /tmp/s11cd-source-trace.XXXXXX)
python -u _measurements/S11c_d_right_source_trace_check.py \
  --left-run /tmp/s11cd-right-source-resume-20260913/left_regression \
  --right-run /tmp/s11cd-right-source-resume-20260913/right_published \
  --run-directory "$s11cdTraceRun/run" \
  > "$s11cdTraceRun/full.out" 2> "$s11cdTraceRun/stderr.txt"
mv "$s11cdTraceRun/full.out" "$s11cdTraceRun/run/full.out"
python _measurements/S11c_d_right_source_trace_check.py \
  --run-directory "$s11cdTraceRun/run" --validate
```

Publication is separate and refuses an existing destination. This rerun uses
the checkpoint's symbolic source operands; it does not run the full d engine.
