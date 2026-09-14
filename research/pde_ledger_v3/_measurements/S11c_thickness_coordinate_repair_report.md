# S11c thickness-coordinate repair: source fix and regeneration

The endpoint-pairing diagnostic is committed at **f796f7ba**. The user then
authorized autonomous continuation and a commit at each substantive step.
This record follows the [repair plan](S11c_thickness_coordinate_repair_plan.md).

The new source trace identifies the discrepancy in **S11c-b's native kinetic
ansatz**. Shared physics sections 1a/1c define the physical thickness through
`e_W = deltaW/W_0` and the kinetic action through its physical time derivative.
The separate local fraction is `(W_0/W_bg)*e_W`, already computed by
`local_thickness_map`. Multiplying that map by background thickness gives
`W_0*e_W`. The native kinetic helper instead uses `W_bg*e_W_t` directly.

The trace extracts the actual density assignment from the native helper,
differentiates it, and compares it with the supplied action after the existing
coordinate map. Both density representatives retain symbolic materials and
background thickness. The extra thickness-inertia coefficient is the computed
`mu_W*(W_bg**2-W_0**2)`. Its retained RIGHT value at the supplied endpoint is
`2*mu_W*W_0**2*eta_bg`. Contracting this difference with the saved independent
frequency legs and row-power maps accounts for the complete saved discrepancy.
No downstream residual has been subtracted or fitted to change a source.

Verification: **142 literal-zero reconstruction/regression scalars**, including
50 entries joining both full five-field inertial pencils to the native source
and 75 entries accounting for all three retained residual matrices. The two
coordinate actions, their nonzero difference, the original nonzero residuals,
and the equal-frequency restriction are separate computed objects. The generic
coordinate identity applies before any numerical witness; endpoint accounting
is explicitly scoped to the saved RIGHT LAB_HELD/RHO4_CONSTANT case.

The [trace transcript](../scripts/out/S11c_thickness_coordinate_trace.out) and
[checkpoint](S11c_thickness_coordinate_trace_checkpoint.json) validate **40
objects and 453 metadata paths**, with restored [L,T,M], explicit epsilon/eta/
sigma and homotopy support, rational numerator/denominator grade metadata for
the coordinate map, and source/packet/output hashes. The output is 146,008 bytes,
SHA256 `932dbe60ef1b6b426a98ee5629bf868f7c7a88c3722706640128c82c91d8557f`.
Measured calculation time after module imports is 21.93 seconds, with 233,392
KiB process peak RSS. The initial metadata-only emission failure is retained in
the scratch `trace` directory; the complete run is `trace-complete` under
`/tmp/s11c-thickness-coordinate-20260914`.

The earlier inertia repair corrected the relative time-variation sign; that
sign remains intact. Its independent action check used the same background-
normalized field expression as the producer, so it did not test this coordinate
choice. The repair will use the existing local-to-physical coordinate map in
the native action and the defining reference-normalized physical field in the
independent test. The frozen historical `committed_strong_rows` stays historical.
No S11b/S11c-a law, authority, face-trace closure or d current formula requires a
change from this evidence.

The source diagnosis is committed at **3c062ae5**. The native kinetic helper now
composes physical thickness through `local_thickness_map` before differentiating
its stationary-background velocity. The existing momentum/time-variation
assembly and the frozen historical comparison are preserved. The independent
time-action check uses the defining reference-normalized physical field and
returns **eight literal-zero inertia-minus-action components** across both
density representatives, with symbolic material/background inputs. The
[action checkpoint](S11c_thickness_coordinate_action_checkpoint.json) pins the
4751-byte [output](../scripts/out/S11c_b_inertia_action_after_coordinate_repair.out),
SHA256 `df08264820d43e3053ef41de07f02eb2d5c5faa35f4340c161492d7cfdfd04db`.
Python compilation and whitespace checks pass.

The source repair is committed at **3b52afcb**. Full b regeneration is running
under `/tmp/s11c-thickness-coordinate-20260914/b_full`; the remaining producers
follow it serially. The serial runner extends the existing
trace-repair staging to b/c1, preserves the established b primaries/single-worker
scope, snapshots inputs, measures resources and leaves publication separate.
Current b/c1/c2/d exports and production outputs still describe the pre-repair
operator until their recorded regeneration completes. Baseline source/export copies
and all accessible S11c output hashes are recorded in
[S11c_thickness_coordinate_baseline.json](S11c_thickness_coordinate_baseline.json).
After the source fix, regenerate b, c1, c2 and d in order, then rebuild actual
endpoint/reference pairing and modal prerequisites. Earlier downstream outputs
remain historical until their replacement is validated. No full S-matrix,
section 3b profile-frequency pole set, global exceptional coverage or d export
is established here. No S10/Lean/authority edit, review/comparator/Wolfram,
downstream physics run or push occurred.

Regeneration checks are prepared for all four b cases (independent action,
baseline source, nonkinetic preservation, full delta accounting and other
slots) and all c2 components (actual imported delta through native field/weak
maps, closure-key independence and other-input preservation). Their numerical
or algebraic results are **not yet available**; Python compilation is the
current check on this preparation. The nine existing d inventories now accept
an explicit plan/report prefix, preserving their historical defaults and data.

The fresh endpoint worker runner rebuilds the physical construction and reuses
only exact operands of the unchanged pure native arithmetic helpers. The cache
reader validated **338 distinct operands and three completed boundary
comparisons** in 8.20 seconds at 78,760 KiB peak RSS. It restores both recorded
stack limits (256 MiB); an initial inspection with an unlimited hard limit was
correctly rejected by the exact runtime guard. This is cache provenance
validation, not a new endpoint balance result. See the
[arithmetic inventory](S11c_thickness_coordinate_arithmetic_cache_inventory.json).
Fresh pairing publication supports a new suffix and never replaces the original
RIGHT discrepancy. End-to-end fresh construction/emission remains to be run
after the repaired production imports exist.

Prepared validation tools are committed at **a01fd555**. A serial continuation
controller now has the concrete b/c1/c2/d command plan and per-stage commit
paths. Its plan-only run and Python compilation complete; actual producer and
physics-check results remain pending. The controller saves its progress under
`/tmp/s11c-thickness-coordinate-20260914/continuation` and writes a tracked
regeneration status at each validated producer checkpoint. It never pushes or
changes a physical source, and it stops with evidence if a validator fails.

The serial controller is committed at **e9b3d6d1** and launched. The
[execution checkpoint](S11c_thickness_coordinate_execution_checkpoint.json)
records the live b producer and waiting continuation, stable current producer
sources and empty stderr at observation. A read-only host process inspection
shows b actively using CPU at 216,720 KiB RSS. The committed baseline b timing
tag records 13,813 seconds (about 3.8 hours); earlier two-hour guidance was an
underestimate for that transcript. The source/export checks and later producers
remain pending. A recurring thread follow-up was rejected by automatic approval
review pending explicit user authorization; none was created. This does not
stop the already-authorized one-time producer/validation/commit queue.

## Export-check execution repair and c2 control propagation

The b native producer completed in **6296.63 seconds**, at **1,942,292 KiB**
peak RSS, with empty stderr and stable source pins. Its saved transcript is
183,361,030 bytes; the new export has 2441 rows and no added or removed keys.
Only `slab_operator` and `slab_operator_term_origins` changed. The production
annex output remains historical until the actual action/export check completes.

The first export check was interrupted after about 98 minutes. Its preserved
traceback locates the cost in `Poly` metadata domain inference, which attempted
a dense GCD over dozens of material coefficients. No physical comparison had
been emitted. The restarted instrument uses the expression coefficient domain
for perturbation polynomials, one projection worker, per-object progress and
timed stack traces. This changes metadata computation, not the action or rows.
The old attempt and traceback remain under `continuation` and `b_checks`; the
fresh attempt is `b_checks_exdomain` under the same pinned run root.

A serialized AST census finds four changed KINETIC provenance components and
28 identical other origin/metadata components. The new algebraic check restores
each actual KINETIC operand and compares it with the independent action and
saved native helper, alongside complete slab-row delta accounting. All other
export values and provenance components have explicit preservation checks.
These checks are running; their final aggregate is not yet available.

The c2 independent conservative-power control repeats the old `W_bg*e_W`
coordinate in its kinetic action. Its source now differentiates `W_0*e_W`,
using the same supplied reference normalization. The c2 full export validation
will emit both source-extracted actions, their differences, current native
assembly joins and imported-origin normalization residuals for all four cases.
It permits a changed KINETIC provenance component only when the completed b
proof accounts for all four of its components; every other provenance component
(including stored energy and face work) must remain identical. These new c2
checks have compiled but have not yet run. The existing face-trace closure is
preserved. No physical authority was changed.

The continuation controller supports a fresh attempt directory and explicit
start stage. It verifies published, committed predecessor artifacts before
resuming, preserves previous logs, and writes `active_continuation.json` for
follow-up discovery. Its c1-through-d plan-only run and Python compilation pass.
The earlier automatic-approval rejection is resolved: the user explicitly
approved continued recurring repairs and commits, and the active 30-minute
thread follow-up is `continue-s11c-d-repair-and-build`.
