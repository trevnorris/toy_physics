# A11/A12 decision boundary after the user's scope critique

2026-09-28. Source-only assessment; no new scientific run, pickle restoration,
review submission or physical binding. **Stop the thickness implementation
chain here.** The next decision concerns nonempty radiation/channel coverage,
not another repair at omega=1. Earlier results and failures remain preserved.

The user's criticism is justified on priority: the closed-slice construction
has not established whether the retained leakage checks are affordable. The
reference-kernel route did have an approved bounded cost/scope decision in the
[constructor plan](../directives/S11c_d_FORM_constructor_plan.md), but that
decision did not establish its applicability or total cost on a radiating
domain. Successive approved ingredient stages did not resolve that question.

## What existing evidence establishes

The [amendment, sections 2 and 3.3](../directives/S11c_d_SCATTERING_FORM_AMENDMENT.md)
distinguishes two obligations: A11 is a thickness-like **end-channel** current
check; A12 is **bulk-depth** outgoing flux. An open acoustic branch alone does
not establish A11, and an open end channel alone does not establish A12.

For the recorded c_s0=10 and tangential components (1/5,1/10), the saved bulk
dispersion gives

    q_depth^2 = omega^2/100 - 1/20 - k_n^2.

Thus nonempty real bulk propagation is algebraically possible when
`abs(omega)>sqrt(5)`, with
`abs(k_n)<sqrt(omega^2/100-1/20)`. At omega=1 it is closed everywhere on the
real normal-momentum line. The threshold agrees with the existing
[frequency-source report](S11c_d_frequency_source_report.md). This is a
kinematic candidate domain, not a selected new input or a validated witness.

The same report's end candidates near 0.27386 and 0.27250 explicitly retain
denominator/opposite-sheet artifacts and are not physical thickness-channel
thresholds. No saved result inspected here establishes the necessary nonempty
flux-carrying A11 branch or a complete A12 witness with incident normalization.
This is not a proof of their absence.

Admissibility must retain the actual constitutive/closure parameters, smallness
and regularity assumptions, and the c1 rest-frame conditions. Away from grazing
these are `abs(q*v/omega)<<1`, `abs(omega*v)/(c_s0^2*abs(q))<<1` and the
independent subsonic condition. Grazing is the stated strict-rest result. The
amendment allows these as conditional validity labels; no invented drain value,
removed physical damping or broad parameter scan is needed to discuss them.

## Reuse and missing work

| Existing work | Transfer to a radiating calculation |
|---|---|
| Frequency-live reduced sources and complete symbolic end pencils | Reusable source inputs. Their actual frequency dependence must be bound; the omega=1 evaluated matrix cannot substitute. |
| Native actions, grade/profile rules, units and face identifications | Reusable structure within the recorded five-field problem. The localized transform implementation has not passed its first actual case, so it is not a validated shortcut. |
| c1 outgoing acoustic fields/flux operands and c2 reference traces/normal jets | Concrete saved starting points, indexed in the [A9 dependency report](S11c_d_A9_dependency_report.md). They do not yet supply the general solved-state-to-bulk current contraction, supported quadratic terms or incident denominator at a radiating input. |
| omega=1 reference inverse, poles, residues and end lifts | Retain as closed-domain evidence. Their evaluated values and input-specific certificates cannot be relabelled for another frequency. |
| omega=1 outgoing prescription | Its Fourier convention and isolated-pole assembly are useful design structure. Its positive radical, nonzero denominator and decaying-bulk current certificates do not cover a real radiating band. |

The last limitation is concrete in
[the outgoing source](S11c_d_outgoing_prescription.py): it uses a positive
`sqrt(a+b*k_n^2)` coordinate with a,b>0, and requires decaying-bulk normalization
for its selected real blocks. Above the acoustic threshold the branch has real
transition points. A radiating action must supply the appropriate boundary
values/continuous contributions there and the corresponding physical flux
normalization. The existing depth-integrated decaying-mode pairing does not
automatically extend to propagating exterior fields. A finite pole sum or a
formal frequency substitution cannot fill this gap.

This need not imply an exterior mesh/PML solver: the existing scalar acoustic
map may permit a direct real-axis Fourier construction. That is a candidate
route to assess, not a supported implementation claim. Keep the complete
coupled rows and actual supported current terms; no scalar proxy or current
deficit may be reported as leakage.

## Bounded cost and recommendation

No measured full radiating-response cost is available. The closed reference
inverse took about 19 seconds, selected historical end-spectrum checks about
19–26 seconds, while an exact closed-slice outgoing continuation took about
an hour. These are different operations and do not predict radiating runtime.
The latter used only about 133 MB, so the record does not support blaming the
2GiB memory cap. Nor has useful parallel speedup of that symbolic operation
been demonstrated. Preserve the containment limits.

The smallest useful next investment is a **single focused admissibility/reuse
decision**, covering the saved frequency-live end relation, actual branch,
physical closure/current and c1/c2 bulk map. Proposed budget: one concrete
method/instrument review packet and, only if source inspection needs actual
saved-expression evaluation, one normally guarded job capped at 900 seconds.
Return a supported candidate domain plus missing objects and an implementation
cost, or an explicit unresolved dependency. No response pilot, sweep or retry
is included in that budget. This is a recommended ceiling, not a forecast that
the physics must resolve within 15 minutes; authoring/review time is additional.

If that assessment shows the existing acoustic map suffices, the later minimal
radiating test still needs three pieces: a supported frequency-bound coupled
outgoing action, its actual face-field/bulk-flux/incident-current contraction,
and one claim-relevant comparison. It is not presently defensible to price
those as a quick rerun. If a new boundary method or unsupported order is needed,
stop with that specific cost/scope choice instead of building it automatically.

My recommendation is to buy only the bounded decision above, not more omega=1
construction. Alternatively, the user can choose a reduced handoff stating
“closed-channel results; radiating leakage not established.” That would be an
explicit scope change; it is not completion of the currently retained
A9/A11/A12 and does not provide a physical leakage magnitude for S11c-e.

## Review and failure status

[G4](../../../CLAUDE.md) says to stop on substantive clearance, not merely two
green labels, but also says a repaired result needs its own review. My local
dispositions are not non-author reviews. Candidate labels and later permission
to run do not themselves discharge that repaired-build acceptance duty. The
historical local-correction/saved-validation evidence remains intact; do not
promote it as fresh independently cleared physics. Any future acceptance must
resolve the applicable repaired-build review or an explicit scoped exception.
This assessment authorizes no review rerun or new submission.

The just-resumed thickness stage stopped at `block17-grade01-local-0-0` with
`unrecognized localized source factor`: 3.691198 worker seconds, 68,558,848 peak
whole-job bytes, zero swap/events, all 23 inputs unchanged. Six restore returns,
eight reported source joins and the incomplete input survive in 37 files.
There is no transformed-source/action result. The failed polynomial guard's
intermediate coefficients were not emitted, so its exact expression-level
cause is not established by the traceback alone. No repair or retry follows.
See the [failure checkpoint](S11c_d_localized_thickness_failure_checkpoint.json).
The completion hook has finished and its event has already been handled.

Future bookkeeping should use concise decision records and existing operation
receipts, without duplicating ancestry inventories solely for reassurance.
Keep necessary input identities, emit-before-guard failures and actual resource
logs. Scratch remains runtime storage and is never committed. Lean and the
shared guard are untouched. No worker is running for this stage.
