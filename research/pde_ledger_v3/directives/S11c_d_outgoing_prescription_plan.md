# Next stage: outgoing reference prescription from saved operands

Status: **method design in progress; no new scientific run or external review**.
The user asked to move on after the accepted reference-inverse stage on
2026-09-27. This continues that work with source and metadata inspection.
The [completed ingredient](../_measurements/S11c_d_reference_kernel_report.md)
is the input, not a reconstruction queue. The scope below is the next bounded
implementation proposal; it does not grant its own method clearance.

## Immediate objective

Specify the missing outgoing real-normal-momentum prescription for the actual
saved full five-field reference inverse. Retain the explicit Fourier density,
normal measure, time sign, inherited radiation branch and physical units.
Distinguish the regular density, local singular blocks and the complete
outgoing integral. The existing physical point remains unchanged.

The immediate next implementation should join the accepted inverse to the
already saved REFERENCE modal/current and acoustic-branch operands, then
construct the supported local prescription pieces. It must identify any
remaining analytic-continuation premise before labelling an outgoing integral
complete. A name such as `G_out` or a bare principal-value integral does not
complete the missing construction.

## Source findings that change the implementation

1. The accepted REFERENCE modal checkpoint already points to
   `_scratch/s11c/s11c-end-normalization-20260914/reference/modal.pickle`
   (434,483 bytes). Its hash matches the checkpoint. There is no reason to
   recompute that mode census.
2. Its recorded inventory contains 18 lifted candidates and 22 basis
   directions: 14 nullity-one records and four nullity-two records. Two
   physically current-normalized subspaces are recorded, each of nullity two,
   with opposite signed currents. These are metadata facts; this inspection
   did not restore or re-evaluate the native packet. A construction assuming
   all real singularities are simple scalar determinant zeros is inadequate.
3. The original channel input and the current development input share the
   original physical parameters and unit frame; the latter adds background
   coefficient bindings. This is a metadata join, not proof of symbol identity.
   The worker must check the actual reference-symbol, coordinate and unit
   mapping. In particular, the checkpoint's `nativeQ` is not silently treated
   as the stored `PHYSICAL_Q` or as the new kernel's momentum coordinate.
4. `continuum_boundary.construct_end` consumes this modal checkpoint. Its
   real-mode selection uses complete signed-current blocks; complex modes use
   outward normal decay, subject to sheet and bulk-decay checks. Those rules
   are useful evidence, but the boundary constructor also builds new mode
   variations and currents and must not be called again for this task.
5. The c1 radiation rule retains outward energy flow or outward decay and
   forbids branch re-selection during frequency continuation. The d current
   convention requires the full left/right normalization and current matrices
   for degenerate spaces. These source requirements take precedence over a
   review suggestion phrased in terms of simple zeros.
6. `EdgeReduction` already saves branch equations, Fourier normalization and
   separate half-line/profile Abel operands. Those operands are reusable.
   The profile regulator is not a limiting-absorption parameter.

The [source receipt](../_measurements/S11c_d_outgoing_prescription_sources.json)
pins the inspected definitions, checkpoints, original/current physical inputs
and exact candidate artifacts. Native value/identity checks remain work for
the guarded implementation; the receipt supplies no scientific clearance.

## Proposed bounded construction

Use only the accepted inverse/context and the own REFERENCE modal packet plus
its explicitly required saved symbolic pairing operands. Import no scientific
producer and run no mode, reduction, current or LU constructor already done.

- Join physical field ordering, normal coordinate, acoustic coordinate,
  frequency, unit frame and full strong symbol to the accepted inverse. Keep
  the actual join operands and residuals. A shared case label is insufficient.
- Recover the inherited branch relation and its selected real-frequency
  boundary value before proposing complex continuation. Do not substitute a
  complex frequency into real-axis Piecewise tests or choose a new square-root
  sign at each point.
- Reuse every relevant saved real-axis mode subspace and its full left/right
  basis, current and normalization matrices. Inspect the already available
  total-normal-derivative operands before constructing any genuinely missing
  local quantity. For a regular semisimple block, test the projected derivative
  matrix before using a block residue; nullity two does not imply a double
  resolvent pole. A singular projected derivative is a specific stop, not
  permission to start an exceptional-point campaign.
- Derive the local outgoing side from the inherited Fourier/radiation/current
  convention and the actual block data. If the full coupled continuation or
  a block's direction is not determined, retain its operands and name the
  missing premise. Do not impose a new global pole-free theorem or infer one
  from the finite mode list.
- Assemble the actual regular and supported singular pieces with the saved
  Fourier measure. Retain the continuous contribution and report explicitly
  whether this supplies an outgoing integral on a stated domain or only
  local prescription ingredients at the existing input.

Required implementation checks are the symbol/coordinate joins, full-subspace
normalization and projected derivative rank, source-derived sign/branch joins,
and a one-sided change to the relevant branch or direction operand that moves
the prescription check. These checks must use actual operands, preserve
returns before guards, and avoid expected physical coefficients.

## Cost, stopping point and downstream work

Proposed ceiling: one new guarded worker, at most 900 seconds, 2 GiB, zero
swap, one CPU, nice 15, 32 tasks and one native thread, using the existing
supervisor and silent local completion hook. This is a cap, not an estimate;
the previous inverse's 19-second cost does not establish the cost here.
The worker, exact consumed inputs and selected checks still need preparation
and the scoped implementation/method gate before a launch. No external review
submission is included in the present source-inspection work.

Stop with saved results on an actual source mismatch, unsupported block,
missing continuation premise or resource cap. No automatic retry, frequency
pole search, new mode census, changed physical input or radiating witness is
part of this proposal. Do not label partial local data a full Green operator.

After a supported outgoing prescription is available, the next substantive
work is the two-asymptote response: unequal profile tails, incident lifting,
forcing/extraction maps, independent eta/sigma grades and actual current
normalization. A11/A12 and all-case FORM completion remain separate retained
obligations. General centre mechanics and the deferred pole queue stay outside
this continuation.
