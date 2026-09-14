# S11c-c2 reference/physical pressure trace repair

User authorized continuation after checkpoint `34a0b38b` isolated the pressure
trace/reference-slot mismatch. Preserve b/c1 inputs and the unrepaired c2/d
baseline under `/tmp/s11c-trace-repair-20260913/baseline` before native edits.

1. Read the actual pressure member of inherited `face_shift`; remove its one
   wave-amplitude factor and extract its reference-value/normal-jet coefficients.
   Derive the multiplication kernel from the existing computed profile Fourier
   definition. Keep the physical trace source and its reconstruction residual.
2. Construct the reference outgoing continuation ansatz. Apply the inherited
   trace to its value/normal derivative to obtain the trace operator. Solve this
   operator in the existing two/three-leg triangular kernel algebra, retaining
   the full eta/sigma rectangle. The normal derivative acts on the output leg;
   retain the mixed composition of the height map with the first-shape response.
3. Differentiate the resulting reference continuation to compute the reference
   normal jet. Solve the inherited affine trace equation for its reference
   pressure slot using that jet and c1's physical-face pressure. Bind these two
   computed reference quantities before closing the full slab operator. Retain
   c1's physical-face pressure and the spectral reference-pressure operand.
4. Emit exact kernel inverse/forward compositions, actual inherited trace and
   source-row pressure/jet factorization operands. A matrix inverse composition
   is an arithmetic check, not an independent physics derivation. Test all four
   anchoring/density cases and both faces, with the mixed grade explicitly live.
5. Build one complete c2 case and its source-work checks, then regenerate the
   four-case c2 export/transcript with the existing semantic export guards.
   Rebuild d from the regenerated inputs, rerun both-end source joins and the
   affected spectrum/frequency/current packets, and publish only complete runs.

No source result is fitted or subtracted to force a zero. Every new residual
is emitted with operands, dimensions and grades before a guard. Retain all
normal-momentum branches and denominators. No new physics premise, S10/Lean
edit, review/comparator/Wolfram/downstream engine is part of this repair.
Only one heavy CAS process runs at a time; preserve partial attempts separately.

## Regeneration sequence after the source joins

All five repair steps are complete, including both source joins, LEFT baseline
regression, full d regeneration and its nine inventories, both end-frequency
refreshes and reference-current source validation/publication. See the
[repair report](S11c_c2_trace_repair_report.md). The reproducible sequence below
records this completed checkpoint. The full native d run is under
`/tmp/s11c-trace-repair-20260913/d_full`. Its producer manifest is written
before launch and completed only after process exit and source/artifact hashes.
Do not rerun into that existing directory or publish a partial transcript.

After the producer completes, run serially from the ledger directory:

```bash
python -u _measurements/S11c_c2_trace_repair_d_recheck.py
python _measurements/S11c_c2_trace_repair_d_publication.py
```

Inspect the actual root/subspace census, sheet-path and exceptional records,
numerical residual/refinement inventories, dimension closure and lossless
codec/index results. Then use the publication instrument's `--publish` flag.
It atomically replaces the main d output and records the preserved annex payload.

For each end, use a fresh directory with `S11c_d_end_frequency_check.py`,
passing the completed d manifest, the existing channel-input file and `--end`.
Validate using `S11c_d_end_frequency_validate.py --publish
--publication-suffix trace_repair`. Rebuild the reference-current source with
`S11c_c2_trace_repair_run_current.py --run-directory <fresh-dir> --manifest
<completed-d-manifest>`, then run the existing nonlocal-current inventory.
Source/payload pins, all computed residuals and metadata must be inspected before
fresh publication. `S11c_c2_trace_repair_current_publish.py --manifest
<current-manifest> --checks <current-inventory>` records the current validation
and publishes to a fresh path. Historical frequency/current packets retain
their old pins.

The both-end physical-current extension is specified in
[S11c_d_both_end_current_plan.md](S11c_d_both_end_current_plan.md). It follows
this validated repair checkpoint; the native producer is no longer running.
