# Saved current balance and numeric-first premise assessment

2026-09-29. Source, existing JSON and byte-hash inspection only. No scientific
payload restoration, numerical contraction, symbolic calculation, implementation,
external review or job launch. The user requested answers to two questions and
at most one proposed job, followed by a go/no-go stop.

## 1. What the saved current records can answer

The four-case finite profile results already contain incident and outgoing open
currents at the original physical contrasts. The finite constructor saves
`incomingFlux`, `outgoingFlux`, `outgoingFluxRatio`, the scattering matrix,
outgoing current matrix and channel labels in `finite-solution.pickle`
([source](S11c_d_finite_scattering.py:387)). The case response also saves
common-origin current ratios in `finite/observable.pickle`
([source](S11c_d_remaining_case_response.py:113)). Saved current selection records
classify the four open directions as transverse at this input; these are existing
selected-channel results, not the still-unfinished new premise test.

The relevant distinctions are:

| Saved object | Present? | What it supplies |
|---|---|---|
| Full finite profile response at the original contrast | Yes, all four cases | Actual incident and outgoing open-current values, including both reflected and transmitted outputs in each column |
| Twelve uniform-background controls | Yes | Same-background homogeneous mode matching and full end-current checks; **not** the same Abel-regulated finite-profile quadrature experiment |
| Zero-contrast coefficient baseline in the regulated continuum system | Yes | `openFractionCoefficients[0]` and corresponding current operands; a baseline of that coefficient system, not a measured additive regulator-absorption term |
| Three contrast-halving field solutions | Yes, all four cases | Direct solutions of the saved coefficient-polynomial systems and their retained expansions; not three new full-unsplit physical profile solves |
| Physical regulator-only loss floor or its subtraction at each contrast | Not established by these saved controls | Neither a physical dissipation fraction nor an error bound is supplied |

The three contrast entries store `direct`, `retained`, `difference`, a maximum
and equation residual. Their producer does not save a fresh full current balance
for each direct solution ([baseline source](S11c_d_continuum_response.py:265),
[other cases](S11c_d_remaining_case_response.py:178)). Existing current series can
support a separately scoped numerical contraction in that same truncated system.
That would be a regulated/truncated current diagnostic, not automatically a
finite-contrast physical total-loss estimate. The saved open-current series
explicitly omit parent pure-second-order terms that can interfere with the
nonzero baseline ([scope](S11c_d_continuum_response.py:245)).

The dedicated uniform controls instead solve a 10-by-10 homogeneous matching
system from saved modes, without the profile quadrature/Abel regulator
([source](S11c_d_uniform_response.py:102)). The four-case report explicitly
distinguishes constant-background distributional symbols from positive-regulator
finite profile quadratures. Their nearly identity transmission cannot be
subtracted as a measurement of the regulator-only floor of the profile runs.

There is a useful **same-system zero-grade baseline**: the existing baseline
response checkpoint records open fractions
`[1.0000000000001514, 1.0000000000005675, 1.0000000000002793,
0.9999999999999827]`. This shows a nearly unit zero-grade result in that stored
coefficient calculation. Subtracting it still does not identify regulator bias
on the changed profile field, or repair the omitted orders/boundary uncertainty.

The regulator is `s11cdAbelRegulator`: it enters the half-line Fourier operands
as `exp(-a*abs(xi))`, before transformation. It is distinct from the physical
permeability and memory closure. Its finite value is not, by definition, an
added dissipative port with a separately known power loss
([native source](../scripts/S11c_d_mixing_scattering_sympy_audit.py:1624)). Thus
“0.1 absorbs this fraction” and an additive regulator-floor subtraction are not
established by inspecting that parameter.

We already have a selected regulator sensitivity measurement; no new reader is
needed to discover it. The saved baseline domain checkpoint reports:

| Existing comparison | Maximum total-current ratio change |
|---|---:|
| 97 to 129 coefficients | `1.8493051801016236e-7` |
| Boundary/source interval 48 to 64 | `2.611544669406385e-7` |
| Abel regulator 0.2 to 0.1 at interval 64 | `1.0949463558063144e-11` |

At the final baseline settings the saved outgoing fractions are
`[0.9999999879053679, 0.9999999307820313, 1.0000000426098308,
1.0000001062580453]`. The four-case report gives an overall approximate range
`0.9999995945` to `1.0000001063`, including the baseline. These values straddle
one and lie within the declared absolute current reporting resolution `1e-6`.
They do **not** resolve a positive physical loss or prove zero loss.

The small regulator change is empirical sensitivity between two finite values,
not the absolute regulator error, a limit, or a bound. Boundary/discretization
effects and retained-order interpretation remain. An apparent squared-contrast
trend would not separate those effects on its own. The fourfold decrease of the
stored field differences is a truncation diagnostic, not a loss measurement.

**Decision on question 1:** the raw current diagnostic is already available at
the original contrast. A new reader could expose more finite-system quantities,
but cannot supply the missing physical attribution or resolution from these
records alone. Do not spend the remaining execution slot on a reader advertised
as extracting resolved total loss by uniform-floor subtraction.

## 2. Binding the fixed point before current-grade extraction

There is a concrete source-supported reordering. The stopped line is
`sp.Poly(x, epsilon_shape).nth(2)`, extracting the quadratic **wave-amplitude**
carrier from the current matrices. It is not extraction of eta/sigma orders
([helper](S11c_d_transverse_face_premise.py:276)). `bound_end` calls this helper
on the unbound current source and only binds its result afterward
([call site](S11c_d_transverse_face_premise.py:408)). Saved inventories confirm
that these matrices still contain material parameters, both frequency legs,
tangential components and normal/depth momentum legs. The RIGHT current matrices
also contain eta; sigma is absent from these particular uniform-end matrices.

For a diagnostic narrowed to `omega=1`, tangential components `(1/5,1/10)`, bind
those exact rational values and fixed material parameters **before** extracting
`epsilon_shape^2`. Keep epsilon explicit for that extraction; retain eta/sigma
independently wherever present until the intended end/grade specialization.
Keep the normal/depth momentum and harmonic-leg distinctions needed by the
existing branch/current checks. Do not set a varying grade to zero or invent an
absent sigma dependence. Save specialized operands and reconstruction residuals.
The substitution map is independent of epsilon, so the proposed ordering does
not require a changed physical model; runtime source/branch checks still apply.

This needs an explicit carried-symbol allowlist: the current `bind` helper does
not accept a raw epsilon carrier unless it is supplied in the extra mapping
([binder](S11c_d_transverse_face_premise.py:257)). Simply moving its call would
therefore introduce another failure. The current frequency/tangent-live guards
also describe a two-dimensional map. A fixed-point instrument must replace those
claims with exact specialization checks and report one-point coverage, not claim
the original 576-record map completed. The earlier map was deliberately symbolic;
binding its fixed material/current context before heavy extraction was still an
available simplification.

**Decision on question 2:** yes, the saved source permits this narrower
numeric-first approach. Its speedup and required runtime cannot be measured from
source alone. Propose a **180-second local cap per end binding** instead of 35,
at most 540 seconds across REFERENCE/LEFT/RIGHT, within the ordinary 840-second
native / 900-second outer envelope. The remaining native time is for the one-point
transverse/current check, limited reference face test, controls and publication.
These are caps, not a promise of completion. Preserve a timeout as unresolved;
no retry or automatic time extension.

## One proposed job, serving question 2

After go/no-go and applicable fresh independent build review, one fixed-point
transverse/current plus reference-face premise job for LAB_HELD/RHO4_CONSTANT at
the unchanged saved azimuth/input. Reuse the eleven completed premise returns;
specialize only the unfinished binding path and existing point/face checks.
No upstream constructor replay, broad map, thickness classification, profile
solve, new loss calculation or additional reader is proposed. Save raw inputs,
completed returns and controls as each finishes. Return established/failed/
unresolved premise statuses and actual cost, then stop.

Keep the shared guard/supervisor, 2 GiB, zero swap, one CPU, nice 15, 32 tasks,
one native thread, one worker and hook-first launch. This would consume scientific
execution **4 of the ceiling of 4** if later approved. Only the earlier executions
were authorized; this note creates no gate or launch authority. A local author
disposition is not independent clearance. No larger method or optional review
campaign is proposed.

This job may establish a pointwise premise; it would still not supply total loss,
all-grade matched-end localization, a completed outgoing field or A11/A12. The
analog-light frequency calibration remains OPEN. **Stop for the user's go/no-go.**

## Inspection provenance

The source checks above were matched to the accepted case checkpoints where
listed. Byte hashes verified 24 selected saved packets: four each of finite
solutions, finite observables, continuum responses, formal remainders, open-current
packets and uniform responses, plus all saved packet sets of the three baseline
domain comparisons. Checkpoint/runtime checks hashes also agree. No pickle was
restored. Exact original routes and hashes remain in
`S11c_d_remaining_case_response_checkpoint.json`,
`S11c_d_remaining_case_flux_checkpoint.json`,
`S11c_d_remaining_case_uniform_checkpoint.json` and
`S11c_d_finite_scattering_domain_checkpoint.json`.
No accepted source, runtime artifact, old failure, Lean file, shared guard or
protected builder suffix was changed by this assessment.
