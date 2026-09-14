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

The export-check preparation is committed at **f618178a**. The restarted b
comparison completed in **595.97 seconds**, at **1,717,600 KiB** peak RSS:
**2634 zero residual scalars**, with **91 objects and 2876 metadata paths**.
The validator restored every emitted record, recomputed carrier fingerprints
and the residual census, checked metadata paths and source snapshots, and found
no nonfinite objects. The four KINETIC provenance changes are fully accounted
for; all 28 other provenance components and other export values are preserved.
The physical comparison remains generic in the retained material/background
symbols across all four cases; this is not a generic spectral witness.

The focused [b export transcript](../scripts/out/S11c_thickness_coordinate_b_export.out)
is 1,361,873 bytes, SHA256
`7620bb53be733b28bfbc1abbd71e940ba1be9ea6db7bde5b436df5a47a4cadd8`.
The complete native b output has been atomically published at its production
path (183,361,030 bytes), preserving its previous annex payload. See the
[b artifact inventory](S11c_thickness_coordinate_b_stage_inventory.json) and
[validated action/export checkpoint](S11c_thickness_coordinate_b_export_checkpoint.json).
c1/c2/d regeneration and fresh endpoint/reference balance remain next.

The validated b artifact checkpoint is **1175b686**. Its native and focused
outputs are both confirmed annex pointers with unchanged payload hashes.
Before starting the remaining queue, the shared export validator now records
improper-integral endpoint nodes separately from nonfinite value nodes. An
actual saved c2 coupling integral has six infinite endpoints and no nonfinite
coefficient; explicit infinite/NaN/complex-infinite coefficient probes and an
infinite-integrand probe remain detected. See the
[domain-syntax checkpoint](S11c_thickness_coordinate_integral_domain_checkpoint.json).
This is validation of expression syntax, not convergence or integral evaluation.
The regression took 72.70 seconds at 3,875,140 KiB peak RSS while selecting the
actual operand from the serialized c2 export. It ran after b validation ended;
no concurrent heavy CAS job was launched. c1/c2/d are ready to resume serially.

The integral-domain validation preparation is committed at **3abb2936**.
The controller resumed at c1 under `continuation_after_b`, after checking the
published and committed b predecessors; its active pointer records the new
state directory. Initial controller and c1 stderr are empty. The older queue
remains preserved as interrupted-attempt evidence.

While c1 runs, the fresh endpoint command plan now lists the actual source,
frequency and independent-frequency worker/emitter/validator invocations, fresh
publication paths and per-stage commit boundaries. The reference-current
publisher accepts a new suffix/checkpoint path and refuses existing targets;
its runner describes the supplied producer rather than the historical repair.
Python compilation, CLI parsing and instrument-path checks pass. No fresh
endpoint source, frequency or pairing result has yet been computed after this
coordinate repair; the native queue remains the prerequisite.

## c2 export-guard restart

c1 completed in **515.16 seconds**, at **1,652,564 KiB** peak RSS. All 44
exported values are unchanged, source/artifact checks pass, and its 90,722,854-
byte transcript is annexed in checkpoint **847d47c6**. c2 then completed in
**1874.03 seconds**, at **2,649,620 KiB**, with empty stderr, stable sources
and a saved 530,883,300-byte native transcript. Its artifact inventory passes;
only `s11cc2ClosedSlabOperator` changes among its 70 exported values.

The focused c2 check stopped at its dependency guard. The emitted b/c1
preservation data contain **2483 zero scalar differences**, but the guard
iterated dictionary keys, so nonempty row names triggered it. The guard now
uses the nonzero count computed from the emitted residual payload. The actual
zero census is accepted and an explicit one-nonzero-count probe is rejected;
see [guard checkpoint](S11c_thickness_coordinate_c2_guard_checkpoint.json).
No native physics source or computed operand changed in this guard repair.

Resume only validation in fresh `c2_checks_values` and
`continuation_after_c2_native` directories. The controller's explicit completed-
producer reuse checks successful exit, stable/current sources, every saved
artifact and the actual current export before skipping the native run. The
stage inventory still runs, followed by physical export checks, publication
and a commit before d. Its plan-only run validated the completed c2 artifacts;
Python compilation and whitespace checks pass. Physical c2 action/control and
closure-delta results are still pending, and production c2/d outputs remain
historical until validated publication. Preserve the failed guard attempt.

The dependency-guard restart is committed at **a866a6e1**. Its fresh check
passed the repaired guard and emitted all four kinetic-control comparisons.
Before publication, inspection found incomplete order metadata on the newly
exposed raw kinetic actions: the inherited structural rule multiplies the
support of a power's base and can omit intermediate orders in a squared sum.
This affects the raw diagnostic action; it does not change the computed action
or the retained native producer. The attempt was interrupted with its output
and traceback preserved in `c2_checks_values` and `continuation_after_c2_native`.

The focused checker now computes exact expression-domain polynomial support
for every small kinetic diagnostic, including the lambda homotopy and zero
operands. Reconstructing **108 actual emitted scalars** under the imported
profile definitions gives **216 zero polynomial reconstruction residuals**;
**50 metadata descriptors** change. See the
[kinetic grade checkpoint](S11c_thickness_coordinate_c2_kinetic_grade_checkpoint.json).
The other closure-object metadata retains its explicitly named native
structural convention. No native source, export or physical operand changed.
The metadata regression took 2.59 seconds at 334,668 KiB peak RSS.

All four native power-residual payloads are byte-identical to the committed
baseline, including their metadata; the
[power-preservation checkpoint](S11c_thickness_coordinate_c2_power_preservation_checkpoint.json)
pins both transcripts and the four exact comparisons. This is preservation of
the recorded raw residuals, not a new simplification or physical verdict.
Resume the completed c2 producer's validation in `c2_checks_exact` under
`continuation_after_c2_metadata`, then let the existing publication/commit/d
queue proceed. The full closure-delta comparison remains to finish.

## Completed c2 validation and annex-content recovery

The exact-metadata restart completed in **364.40 seconds**, at **2,632,892 KiB**
peak RSS. All **2846 residual scalars are zero**, across **418 objects and
3244 metadata paths**. The 44 component records account for all four cases;
only `s11cc2ClosedSlabOperator` changes among the 70 native export values.
The native output, focused transcript and export/inventories were committed at
**537d78fd**. See the
[c2 export checkpoint](S11c_thickness_coordinate_c2_export_checkpoint.json).

The controller's post-save hash check then stopped the queue before d. The
native annex object contained 79,560,704 bytes, an exact prefix of its expected
530,883,300 bytes. The complete validated producer output remained intact in
`/tmp/s11c-thickness-coordinate-20260914/c2_full/full.out`. Disk space was not
exhausted; the cause of truncation has not been established.

Targeted `git annex fsck` quarantined the bad object. `git annex reinject`
restored the existing key from a verified recovery copy, preserving the original
producer output and the quarantined payload. No Git pointer, export or physics
source changed. SHA256 and file-size comparisons now match all three native
b/c1/c2 outputs and both focused transcripts; targeted annex fsck succeeds for
all five. The recovered native c2 SHA256 is
`9712191e3af5e7bbc2eb65824ca283f0c5cb3a35d5e21632d0821412ec018864`.
The [recovery checkpoint](S11c_thickness_coordinate_c2_annex_recovery_checkpoint.json)
records the observed failure, quarantine, reinjection and verification results.

Resume the existing serial controller at d in a fresh
`continuation_after_c2_annex` directory after committing this recovery record.
Its predecessor and post-save content checks remain enabled. Fresh endpoint
sources, independent-frequency pairing and full current/adjoint normalization
remain pending; c2's action/export checks alone do not settle that balance.

Recovery is committed at **5a81223f**. The controller resumed at d at
**2026-09-14 15:37:45 UTC**, after verifying its committed predecessors.
The scoped-symbol worker was observed live using CPU, with empty controller
and worker stderr. Full d production, rechecks, source joins and publication
follow serially; no fresh endpoint check has started.

## Native regeneration complete; fresh endpoint construction started

The full four-case d producer completed in **4166.19 seconds**, at
**2,476,192 KiB** peak RSS, with stable source pins and empty stderr. All nine
recorded inventories, scoped/full source joins and publication checks passed.
The production checkpoint is **8b2e3cf2**. Its annexed transcript contains
**83,989,687 bytes**, SHA256
`19822837896b084fb8bcbd11b16aa5120e2c5b01c1c3084357aed279b7b70ce8`;
the controller's post-save hash check and a fresh targeted annex fsck pass.

The native spectrum inventory retains 24 packets and 432 root/lift candidates,
including 336 scalar and 96 two-dimensional nullspaces. Joint-sheet records
contain 480 transported paths and 24 separate branch-locus-on-path records.
The numerical residuals, exceptional and threshold records remain in their
[full inventory](S11c_thickness_coordinate_d_full_checks.json). These measured
domains do not establish global coverage or profile-frequency bound poles.
All ten broad native outstanding constructions retain their current scope.

The prepared endpoint plan is now executing. Fresh RIGHT source construction
started at **2026-09-14 17:14:12 UTC** in `end_source_right`, using the new full
d manifest and unchanged supplied input. The active endpoint operation is
recorded at `/tmp/s11c-thickness-coordinate-20260914/active_endpoint.json`;
the native `active_continuation.json` now identifies a completed queue, not a
live job. Follow the endpoint operation before launching any other heavy CAS
job. Validate, publish and commit its completed packet before the next planned
stage. LEFT/reference sources, both frequency packets and the three fresh
pairings remain pending. The original retained two-frequency discrepancy is
not yet retested by these successful native regeneration checks.

RIGHT source construction completed in **453.55 seconds**, at **201,604 KiB**
peak RSS. Its saved census has **127 cancellation identities**, all zero, and
no retained nonzero source residual. Its validator/publication remains next;
these construction results do not yet supply a fresh two-frequency pairing.

The prepared endpoint stages now have a serial executor,
`S11c_thickness_coordinate_endpoint_continue.py`. It can adopt the completed
RIGHT construction without repeating it, runs each existing validator, checks
published payload hashes and the recorded residual census, and commits before
the next stage. It uses an exclusive controller lock, source/input/plan pins,
fresh log files, and stops on any command or publication failure. A failed
pairing emitter preserves its full diagnostic output before the stop. The
separate `active_endpoint_controller.json` points to its live progress; consult
it before the earlier single-construction `active_endpoint.json`.

Compilation and the complete eight-stage command-plan check pass. Isolated
operational fixtures exercise accepted exact payloads and rejection of nonzero
residual counts, missing census fields and changed payloads; see the
[controller checkpoint](S11c_thickness_coordinate_endpoint_controller_checkpoint.json).
These are execution guards, not physical evidence. No native construction or
validator changed in this preparation.

The executor preparation is committed at **195566c2**. The serial endpoint
queue started at **2026-09-14 17:24:43 UTC**, adopted the completed RIGHT
construction and began its validator/publisher. Runtime discovery is
`active_endpoint_controller.json` under the run root; the tracked execution
checkpoint explicitly gives this pointer priority over the two completed
earlier queues. Subsequent stages retain the existing order and separate
publication/commit boundaries.

RIGHT source validation and publication completed at **fba5f1fb**. The
validator checked **85 source objects**, **368 source metadata paths**,
**33 retained source residual scalars** and **3484 tags**; all retained source
residuals and all 127 cancellation identities are zero. The annexed transcript
has **8,788,422 bytes**, SHA256
`3619e8ea9da541706ffe14bc9093037104013a942fb85a2ede275c0dbac84501`.
The post-save hash check succeeded. LEFT source construction started
automatically at **2026-09-14 17:27:24 UTC**.

Read-only inspection located the remaining reference-specific loaders,
zero-background bindings and per-mode emitter labels in the modal/adjoint
instruments. Their required endpoint adapters are recorded in the active
repair plan's normalization handoff. The native constructors already retrieve
the endpoint energy/acoustic context. No live pinned implementation changed;
normalization work still requires the completed fresh pairing checkpoints.
