# Numerical radiating pilot: source-only feasibility

2026-09-29. Assessment requested after the final RIGHT stop. **Conditional GO
for a bounded numerical pilot at omega=3; the existing scripts cannot be run
there unchanged.** No implementation, scientific import, payload restoration,
new calculation, or review submission was performed for this assessment.
This is a planning recommendation, not a reviewed method or a launch authority.

The exact-symbolic omega=1 premise track is parked. Selected REFERENCE/LEFT
support remains at the scope of execution5; RIGHT is unresolved because the
restored-summary guard stopped execution7 before its physics checks. The two
saved sourceFactor strings differ by algebraic presentation. No further
summary diagnostic or exact-track repair is queued.

## What can be reused

The actual finite solver already assembles all five fields, 80 nonlocal rows
and 35 source amplitudes, solves four incident columns, and contracts physical
current matrices including cross terms. Its finite-domain wrapper already
varies interval/source extent and the Abel regulator. The required observable
is therefore close to an existing output, rather than a new exact Green
operator or Born construction.

There is also a completed frequency-dependent version, not just a prospective
one: the accepted frequency-matrix checkpoint records 645 unknowns, all 80
rows, 177,848 new momentum nodes and 88.503 worker seconds at `1-0.01i`.
This proves numerical frequency rebinding and assembly have worked locally;
it does not establish a radiating calculation at 3. The saved full frequency
sources and LEFT/RIGHT pencils exist and their selected hashes match.

The bulk has already been eliminated through its outgoing boundary relation.
The native c1 code gives the flat pressure/velocity factor `rho_m*omega/q`
and an outgoing layer factor `exp(i*q*s)/(i*q)`; c2 binds the positive-real
frequency outgoing/decaying square-root branch on each Fourier leg. A new
three-dimensional exterior mesh is not the first route to try. The numerical
pilot would approximate the Fourier integrals and finite in-slab ends while
retaining that physical bulk relation and the physical permeability/memory.

## Changes that actually matter

| Required change | Concrete source evidence and limitation |
|---|---|
| Bind the original frequency-live sources at each real frequency, with numerical material/tangential inputs first. | Every one of the 80 rows depends on frequency. The omega=1 matrices cannot be reused as omega=3 matrices. The later analytic chart is explicitly restricted to `abs(omega-1)<=1/4`; its `i*sqrt(-radicand)` expressions and nonzero-denominator certificates must not be extrapolated to 3. Use the original real-frequency branch operands. |
| Resolve the acoustic branch endpoints in the numerical integration. | At the saved tangents and sound speed, the recorded relation is `q_depth^2=omega^2/100-1/20-k_n^2`. At 3, the real-depth band is `abs(k_n)<1/5`. Current momentum panels target transfer/Abel concentration, not these endpoints on all input/output/middle legs. Complete coupled coefficients must be checked for endpoint cancellation or the appropriate integrable boundary value; merely avoiding an exactly zero node is insufficient. An unsupported singularity is a stop, not permission to invent a prescription. |
| Supply numerical outgoing end maps at the new real frequencies. | `finite_scattering.boundary_map` rejects candidates without `BULK_DECAY_DISK_CERTIFIED`. The old matching-current constructor also requires positive imaginary depth momentum for its infinite-depth integral. Existing numerical subspace continuation/map assembly is reusable structure, but saved modes, current normalizations and seed labels do not transfer to 3. Do not drop the decay test and call every complex root outgoing. |
| Contract only reflected/transmitted transverse output against transverse incident current. | The old solver sums all modes labelled `open`; that happened to select transverse modes in the saved omega=1 run. At a new frequency, select and verify the transverse subspace explicitly, preserve its full current matrix, and numerically check its end drives/current. No infinite-depth decaying-bulk formula is applied to a propagating bulk wave. |
| Check numerical resolution at the new frequency. | The old 129 coefficients per field resolved omega=1. More oscillatory fields at 3 may need a larger basis. Domain/regulator comparisons alone cannot diagnose unresolved waves or acoustic endpoint quadrature. Momentum cutoff and quadrature/basis resolution belong in the small sensitivity budget. |

The main uncertainty is whether the existing finite modal end approximation
is adequate once the bulk continuum is open. Moving the ends and refining
the numerical integration can expose inadequate behavior, but cannot prove
a transparent nonlocal boundary. If adequate end maps or endpoint treatment
require a new exterior-boundary theory, this stops being the cheap pilot.

## Useful deliverable and bounded proposal

Start with **LAB_HELD/RHO4_CONSTANT, omega=3, saved tangents (1/5,1/10)**
and the saved material/profile/contrast values. A possible small neighborhood
is 2.8 and 3.2; these are proposed sample frequencies, not inspected results.
Do the central point first. Add neighbors only if its boundary, resolution
and regulator sensitivity leave an interpretable signal.

Report, for each selected transverse incident column,

    D_T = 1 - (J_transverse_reflected + J_transverse_transmitted)/J_incident.

This is a finite-model current deficit. It can support a rough total removal
from transverse propagation if the incident/end current and numerical checks
hold. It does not by itself separate bulk radiation, thickness conversion
and face absorption. Reflected transverse light is surviving output. A
positive regulator can contaminate the deficit; a difference between two
regulators is not its absolute effect or a subtractable absorption floor.
Retain signed values and report unresolved if the signal is comparable to
the observed sensitivity. Bulk availability at 3 does not predict large loss.

A practical proposed budget is roughly **6–10 finite solves**, including
the central result, a same-setting uniform control, selected resolution,
domain and regulator comparisons, and at most two neighboring points. The
uniform control diagnoses numerical/background behavior; it is not assumed
to cancel the profile-dependent regulator error. A failed or unsettled
central result ends the pilot rather than triggering a broad sweep.

Estimated effort: **2–4 focused working days for implementation and the two
method/build reviews, then roughly 1–4 hours of numerical work** if branch and
boundary handling fit the existing formulation. These are low-confidence
planning estimates, not measured omega=3 timings. Evidence for the scale is
the saved 49.46/51.33/60.61-second omega=1 domain cases and the 88.50-second
nearby complex-frequency construction. Finer bases and endpoint integration
may substantially increase cost. A proposed hard preparation stop is two
authoring days: if a concrete bounded numerical boundary/integration method
is still missing, return no-go rather than start another symbolic campaign.
Substantive review blockers likewise return to the user; elapsed review time
is not a promise of clearance. None of this budget is launched or approved
by this assessment.

The result would be a numerical toy-model estimate with an empirical
sensitivity envelope, not a certified error bound or full FORM/A11/A12
clearance. The physical analog-light frequency band remains uncalibrated.
No exact RIGHT-premise completion, old candidate Green construction or
Born/full anchor is a prerequisite for evaluating this separate numerical
proposal. Its own numerical end/boundary/current checks remain necessary.

## Existing result and process decision

The saved omega=1 finite calculation resolves no loss at its declared
`1e-6` current resolution. Its four baseline ratios describe incident
columns, not four cases. This is not a physical upper bound. The omitted
pure-second-order interference terms in the retained series must not be
automatically attributed to the full finite original-contrast solve; finite
boundary, regulator, discretization and model-truncation limitations still
apply to the latter.

Under the user's new direction, pure serialization/comparison/formatting
repairs receive ordinary local verification and a short non-physics rationale,
not their own two-reviewer cycle. A repair that changes accepted values,
mathematical tests, equations, radiation treatment or claims is not purely
tooling. New method/build work still gets fresh Claude and Grok review.
This record changes no scientific method or acceptance claim and initiates
no review. Stop here for the user's decision.

## Inspected evidence

Source reads and JSON/hash inspection only; no scientific payload was opened.
All six saved domain-case system/solution files, four selected frequency-source
files and three selected complex-frequency system/solution/end-map files
matched their recorded hashes. Saved worker times above come from checkpoint
JSON, which differs slightly from report-level elapsed timing.

- [Finite assembly and current contraction](S11c_d_finite_scattering.py), `load`, `boundary_map`, `construct`; [domain wrapper](S11c_d_finite_scattering_domain.py), `CASES`, `adapters`, `case`.
- [Frequency sources](S11c_d_frequency_source.py), `sources`, `end_sources`; [source checkpoint](S11c_d_frequency_source_checkpoint.json), [source report](S11c_d_frequency_source_report.md).
- [Frequency chart](S11c_d_frequency_chart.py), `root_chart`; [finite frequency builder](S11c_d_frequency_matrix.py), `bind`, `end_maps`; [matrix checkpoint](S11c_d_frequency_matrix_checkpoint.json).
- [Numerical end maps](S11c_d_frequency_end.py), `maps`, `continue_pair`; [native engine](../scripts/S11c_d_mixing_scattering_sympy_audit.py), `BoundedActionQuadrature.compiled`, `FiniteMomentum.rule`, `ThreeMomentum.batches`, `TwoEndedMatchingChannels.construct`.
- [c1 boundary relation](../scripts/S11c_c1_bulk_closure_sympy_audit.py), `dtn_flat_symbol`, `hanzawa_first_kernel`; [c2 branch binding](../scripts/S11c_c2_selfenergy_fold_sympy_audit.py), `outgoing_spectral`.
- [Finite-domain checkpoint](S11c_d_finite_scattering_domain_checkpoint.json), [saved-current qualifications](S11c_d_transverse_face_fixed_point_saved_current_note.md), [final RIGHT stop](S11c_d_transverse_face_right_final_failure_report.md).
