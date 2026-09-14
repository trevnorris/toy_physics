# S11c thickness-coordinate repair: source diagnosis

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

The producer repair and regeneration are next. Baseline source/export copies
and all accessible S11c output hashes are recorded in
[S11c_thickness_coordinate_baseline.json](S11c_thickness_coordinate_baseline.json).
After the source fix, regenerate b, c1, c2 and d in order, then rebuild actual
endpoint/reference pairing and modal prerequisites. Earlier downstream outputs
remain historical until their replacement is validated. No full S-matrix,
section 3b profile-frequency pole set, global exceptional coverage or d export
is established here. No S10/Lean/authority edit, review/comparator/Wolfram,
downstream physics run or push occurred.
