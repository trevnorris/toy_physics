# Four upstream repairs: delivered review and source disposition

2026-10-01. **Composition not cleared; second independent verdict missing.**
Claude delivered a scoped report. Grok emitted only its opening sentence and
ended with `stopReason: cancelled` after one turn (11.278 seconds), despite CLI
exit 0. The reason for that cancellation is not established by its stderr.
It is neither a CLEAR nor a NEEDS REVISION verdict. No retry or new export has
been made. Both processes finished before this adjudication; no peer report was
shared.

The [record](S11c_upstream_repair_review_record.json) preserves both literal
texts, receipts/stderr, every immutable reviewer probe source and complete
stdout/stderr, resource records and file hashes. All 643 inspection checks
passed for the approved 86-file packet, archive/private copies, source excerpts,
historical diffs, consumed pins and runtime evidence. Packet SHA256 remains
`7cf2a67660fcec7f4029fdef6384a0520104b265596f902b2f3486f5e87eee29`.
This inspection read source, JSON and opaque bytes only; it did not import CAS
or restore scientific payloads.

## What Claude actually established

| Repair | Literal scoped conclusion | Evidence and remaining limit |
| --- | --- | --- |
| Inertia sign | Clear on stated domain; fold ordering not computed | Independent action and four source-extracted native functions agree for gradient-free test energy, stationary symbolic W_bg, constant mu_W and both density rules. Full material-constraint fold remains source-inspected only. |
| External work | Orientation clear; non-flat face-row content insufficient evidence | Native multiplier and independently derived flat thickness force agree with responsive sign mutation. Actual non-flat u/centre routing and c2 power pairing were not independently executed. |
| Pressure trace | Two-leg clear; three-leg mixed grade needs revision | Eight face/case two-leg samples agree, with responsive alternatives. One LAB_HELD/RHO4/Plus three-leg sample exposes a conditional missing-direct-term concern; the other seven three-leg combinations were not supplied. |
| Thickness coordinate | Clear on stated domain | Native kinetic rows use deltaW=W0 eW; reinstating W_bg eW produces the stated symbolic nonconstant-background residual. No new varying-mu_W, advection or centre-kinetic result. |

These are one reviewer's findings, not paired clearance. The two-leg numerical
sample's printed zero is not promoted to a global symbolic identity. The current
near-unity restriction, face/current and grazing checks still have to be made;
none follows from these repair checks.

## Mixed-grade finding: credible, but not yet a benchmark correction

The source part is confirmed. In
`scripts/S11c_c2_selfenergy_fold_sympy_audit.py:400–418`, `z_three[0,2]` is zero.
Consequently the mixed response contains the ordered product
`-a G Z1 G Z1 G`, but no direct `G Z2 G` contribution. `build_face` consumes
that response. `retained_shape` at line 664 retains eta*sigma_W, and c1's
governing multigrade section explicitly keeps first order in each bookkeeper.
Therefore a genuine height-times-slope term cannot be dismissed merely as
"second order." The reported issue concerns completeness of the inherited
curved-face composition; it is not evidence that the ordered trace inverse was
implemented with the wrong sign or order.

The missing applicability join is also real. S11c-a defines
`W_bg=W0[1+eta*w1(y/L_W)]` and `sigma_W=eta*W0/L_W`; changing the independent
bookkeepers through eta and L_W must retain the profile/jet relation. Reviewer
`p3` instead fixes `d(x)=cos(q*x)` and q=0.8 while independently varying the
height eta*d and the slope sigma*d'. Reviewer `p2` assigns independent numerical
functions to Fourier profiles and their jets. Those are useful formal probes,
but their relation to the actual two-scale source family has not been shown.
This does not disprove the finding or forbid formal independent bookkeeping;
it makes the full source/profile-jet join the next mathematical dependency.

The report calls p3 an "exact benchmark." Its implementation is a 15-mode,
96-node, 40-digit finite lattice calculation using one mixed finite-difference
step h=1e-6. Agreement with the reviewer's analytic formula is approximately
3–9e-13 at its nonzero pressure entries. There is no mode/node/step convergence
study, so this is not an exact continuum result or an error bound. At the p2
sample the direct-term discrepancy is about 3.608e-5; it is an operator-sample
difference, not a leakage fraction or an error estimate for our finite balance.

**Disposition: retain as an unresolved, substantive mixed-grade consistency
finding.** Before changing equations, derive the candidate direct term with the
actual two-scale profile and Fourier-jet definitions and join its grade and
source addresses to the native c1/c2 operands. Do not quietly add it as a
tooling fix, declare the benchmark numerically wrong, or regenerate producers.

## Coverage and provenance qualifications

- Native `constraint_fold_from_source` exists in b at line 2494 but its body was
  omitted from this packet. Supplied `build_operator` shows kinetic subtraction
  before the fold and readdition afterwards. Full fold validation is uncomputed,
  not evidence of a missing implementation.
- The non-flat geometry objects also exist locally. S11c-a's producer has face
  maps at line 841, traction at 919 and work construction at 953–987. Its current
  `S11c_a_exports.py` contains `face_velocity` at 3131, `traction` at 14119 and
  `virtual_work_shape_deriv` at 14778. A later bounded source extract can supply
  them without a producer replay. They were not literal operands in this packet.
- The old mechanical checks do pin older b/c2 revisions and the obsolete
  W_bg kinetic coordinate. However, the packet also includes the newer
  `S11c_thickness_coordinate_b_export_checkpoint.json` and corresponding c2
  checkpoint, which pin current producers `666b4005…` / `0e90c751…` and current
  exports `3c6555f9…` / `2ba5ed48…`. Their action/delta/nonkinetic-preservation
  records are prior author evidence. They establish current provenance, not a
  substitute for the missing independent fold/non-flat tests.
- The stiffness anchor fixes the stored-row convention; flipping its supplied
  energy on both routes cannot independently test the external-work sign.
  c2's power construction shares b's supplied energy/constraint data. Preserve
  that dependency rather than calling it an independent energy derivation.

## Probe execution and preservation

Claude ran 11 small probes: seven completed and four failed during probe
development (bad import, unresolved float conversion, symbolic-to-complex
conversion and an overstrict derivative stub). The original scripts, partial
outputs and exceptions remain intact. The successful replacement for the last
stub changes no native body: an unused thickness gradient was being asserted
zero too early. Its corrected stub remains valid only for that gradient-free
test. These are reviewer attempts, not retries of any production calculation.

Total supervisor worker time was 224.527 seconds, maximum sampled cgroup peak
68,513,792 bytes (65.34 MiB). Every probe verified 4 GiB cgroup **and** native
address-space caps, zero swap, one CPU/thread, 32 tasks, 4 GiB host reserve and
the 16 GiB aggregate pool. All sampled memory events were zero; every unit had
`RuntimeMaxUSec=infinity`, `Restart=no`. Successful scientific stderr is empty;
failed stderr is preserved. Three denied direct shell commands are recorded in
the report metadata; the actual probes used the required dispatcher. No new
probe, reviewer invocation or scientific restoration occurred in this inspection.

## Consequences for the two tracks

The mixed shape/gradient concern and missing non-flat routing need not become
prerequisites for the **constant-end selected uniform** question. The proposed
term should disappear under a flat, constant profile specialization; that
disappearance must be checked from actual source operands, alongside the full
transverse lift. The inertia/thickness evidence is relevant support, not a
uniform result or full dependency clearance. Keep the uniform findings
provisional and address its separate NEEDS REVISION method findings.

The nonuniform retained composition and defect benchmark are not cleared. No
new loss interpretation, defect sweep, producer replay, equation change or
external submission follows from this event. The concrete remaining decisions
are the corrected uniform method, a source-consistent mixed-grade check, the
missing non-flat/fold evidence, and a genuinely delivered second independent
assessment. They are not an automatic review or repair campaign. Preserve the
rest-bulk/calibration/flow limits and all prior results and failures.
