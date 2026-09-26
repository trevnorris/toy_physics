# Centre-load compatibility: completed bounded diagnostic

**The necessary face-load compatibility check passes in all four cases.**
The new calculation finds zero centre generalized load when the existing
retained c2 thickness-driven closures are used. This is a useful source-level
result, not a complete elimination or uniqueness proof for independent centre
motion. The bounded task is complete; no response or radiation-method job is
launched by this disposition.

The [checkpoint](S11c_d_centre_compatibility_checkpoint.json) accepts the saved
**diagnostic record only**, with its original operands, actual operations,
controls, source qualifications and final resource evidence. It does not give
independent physical clearance to c2 or to a general centre reduction.

## Computed result and its provenance

Let `Fplus` and `Fminus` denote the actual normalized thickness face-load
contributions consumed by c2, and let `S` and `D` denote its saved half-sum
and half-difference. Each case independently produced:

```text
rplus  = -2/W_0
rminus =  2/W_0
centre non-slot remainder = 0
Fcentre = (rplus+rminus) S + (rplus-rminus) D = (-4/W_0) D = 0
```

These are transcriptions of the saved computed returns, not supplied expected
values. Both pressure and global-normal-jet relations passed on each face.
The actual open b `E_W` pressure slots joined the normalized face row, and the
whole b export was previously joined to the successful c2 consumer. The only
scalar remaining in the relation is source-defined positive `W_0`; its inverse
has length dimension -1 and maps the saved thickness-row units to the saved
centre-row units. Cancelled discovery factors do not restrict parameter
branches: the reduced relation was checked against both complete slot
coefficients. The non-slot thickness remainder was retained separately.

The zero uses the existing literal `D=0` return. It is **inherited from that
closure result**, not a second derivation of the closed force or an independent
check of the old acoustic face signs. The applicable c2/c1 fidelity and
cross-engine debts therefore remain. The reused rows are after closure,
profile replacement, retained-shape projection and physical-field mapping,
before `extract`. Formal integrals and the actual case, density, units,
independent grades and source conventions are preserved. No old closure,
Taylor projection, derivative, integral, mode or response was recomputed.

| Own case | New recorded calls | Exact preceding new-call uses | Centre image |
|---|---:|---:|---|
| LAB_HELD / RHO4_CONSTANT | 132 | 257 | saved computed zero |
| LAB_HELD / RHOBR_CONSTANT | 85 | 49 | saved computed zero |
| MATERIAL_ADVECTED / RHO4_CONSTANT | 85 | 49 | saved computed zero |
| MATERIAL_ADVECTED / RHOBR_CONSTANT | 85 | 49 | saved computed zero |

These are recorder calls, not a count of SymPy's internal primitives.

## Controls and saved inspection

The selected whole-face sign change in the imported centre row passed the
local relation checks and produced a final expression structurally different
from the baseline: its computed sum weight is `4/W_0` and difference weight
is zero. The entire result is saved. Its generic `is_zero` property is
**undecided**; no universal nonzero claim is made.

The one-pressure-sign change failed the safe relation gate because its jet
residual was not established zero (`is_zero=None`). The exchange of the two
jet slots was degenerate on this actual input: the changed row equals the
original. That control supplies no additional discrimination. These are
imported-routing controls, not independent acoustics validation.

Final saved inspection checked all 791 operation receipts and their complete
serialized arguments, context, values and exact reuse owners. It traversed
pickle syntax with opaque constructor tags: no saved callable, scientific
codec, symbolic operation or proof ran. All 43 consumed file routes and 2447
manifested artifact files matched their hashes. The two manifest/check files
bring the complete directory total to 2449. This is bounded saved-data
inspection, not an independent mathematical recomputation or ancestry audit.

The worker took 25.648 s, the native supervisor 26.120 s and the guard 26.246 s.
Actual guard/supervisor/child exits were zero, strict stderr was empty, and
checks/stdout bytes matched. Peak memory was 802721792 bytes; memory-cap,
OOM and OOM-kill events and swap were zero. The actual limits were 900 s,
2 GiB, zero swap, one CPU, nice15, 32 tasks and one native thread. Source,
input, logical/canonical and protected builder-suffix hashes stayed unchanged.

## Precise remaining upstream dependency

This closes the question **whether the recorded thickness-driven face closure
itself leaves a centre-force obstruction**. It does not select the independent
centre state. The consumed b kinetic and constraint constructions carry the
in-plane/thickness fields and eliminate density variation; they do not supply
a complete independent centre equation or its solution selection. A homogeneous
centre response can remain possible even when the thickness-driven forcing
vanishes. No new inertia or missing physical coefficient is assumed here.

Before the five-field state can be advertised as determining both general
physical-face drives, upstream b/c1/c2 must supply one concrete disposition:

- the complete centre balance or constraint, joined to the independent
  `ZETA_C` face direction, together with the prescribed centre/incident drive
  and boundary or initial conditions that select its admissible solution; or
- an explicitly prescribed independent centre drive, with a source-supported
  statement of the resulting restricted scattering problem.

If a homogeneous centre sector is excluded, that exclusion needs its actual
operator/domain and solution-selection justification. The present cancellation
alone does not provide it. It also does not show that a new centre solver is
necessary, or that the existing five-field results are wrong.

This is the stopping point required by the approved task and amendment §2.
Do not impose `zeta_c=0` inside d, add a field, or begin an upstream rebuild on
this record. The separate generic radiating-boundary/continuous-spectrum
method remains unapproved and unresolved. A9/A11/A12 completion, the retained
export and independent/Wolfram/T7 duties remain; Option B pole deferral is
unchanged. The next scope decision should address this specific upstream
centre-state disposition before commissioning that larger method.

## Execution references

Implementation `8470e88704a5d93f8753be7ff2829a7da6e64fe0`; launch
`7a6dd237cc0c27d51ed546ef949061591ac6280d`. Original reviewed plan and both
NEEDS REVISION reports remain at `a0de87f1`; bounded correction closure is in
the [adjudication](S11c_d_centre_plan_adjudication.md), with no invented new
independent CLEAR verdict. Metadata preparation/harness failures are preserved.
The producer checks SHA is
`c3145347bd889203d831c03f4c96a8cbd74766d3d8cc15e44712d6a690ba9099`.
Full routes, guards and saved-inspection evidence are referenced by the
checkpoint; no scientific corpus was copied or rewritten.
