# S11c-d: bounded centre-motion compatibility check

Status: **focused gate corrections and local instrument checks completed; ready for guarded launch**.
The user approved the bounded task in the
[boundary disposition](../_measurements/S11c_d_FORM_boundary_disposition.md).
This plan implements that task under the accepted
[Option B amendment](S11c_d_SCATTERING_FORM_AMENDMENT.md), especially §2's
slab-to-face dependency rule. It does not revise the amendment or authorize
the general radiation method. Codex is the author.

## Question and stopping point

Determine what the existing centre generalized force says about the actual
five-field thickness-source sector consumed by c2/d, separately for
`LAB_HELD`/`MATERIAL_ADVECTED` and `RHO4_CONSTANT`/`RHOBR_CONSTANT`.
The deliverable is a necessary compatibility diagnostic and an explicit
upstream dependency disposition, not a new scattering solution.

S11c-b §1a and c1 §1a retain the independent centre displacement. The saved
b kinetic construction has in-plane and thickness velocities; its saved
constraint fold eliminates the density virtual displacement. Neither fact
by itself eliminates the centre displacement. The c2 source identifies the
face velocity with the saved `DELTA_W` direction. The check must distinguish
these supplied/consumed choices from a derived invariant-sector result.

The source review must account for the supplied energy, constraint and face
work before interpreting a residual. An identically vanishing face load is
not a complete centre dynamics or uniqueness proof. A nonzero load on
arbitrary trial fields is not automatically a nonzero load on the solved
five-field sector. Report these distinctions; do not solve a new response
or infer that all previous results need rebuilding.

Stop with a usable source-supported restriction/reconstruction **only if**
the available premises establish it. Otherwise identify the missing centre
equation, prescribed drive, initial/boundary condition, or solution-space
compatibility test and the claims that depend on it. Do not supply missing
physics, impose centre zero, add a field, or extend the job to resolve it.

## Exact saved inputs and their role

The navigation manifest accompanies this plan. It records logical/canonical
paths, literal record/field byte ranges and hashes. The preparation scripts
use only lexical parsing of existing output and literal key decoding; they
have performed no mathematical validation. Full scientific values remain at
their original addresses.

| Input | Actual source / saved record | Use |
|---|---|---|
| Independent face directions | a `dof_fields`, `build_material_face_source`, `supplied_face_maps`; `PY_S11CA_FACE_VELOCITY` and `PY_S11CA_VIRTUAL_WORK_SHAPE_DERIV` | Distinguish the two directions, outward conventions and each anchoring. The material direction records contain shared background advection; adding them wholesale would double that term. No new velocity construction. |
| Centre work and normalized slab face load | b `face_generalized_force_rows`, `mechanical_work_row_normalization`, `build_operator`; `LOCAL_SLAB_FACE_GENERALIZED_FORCE_ROWS` and `SLAB_OPERATOR_TERM_ORIGINS` | Restore the saved centre row and actual normalized `FACE_VIRTUAL_WORK/ROWS/E_W/EXPANDED`. Do not replace the latter with an unnormalized virtual-work coefficient. |
| Kinetic and constraint context | b `kinetic_balance_from_energy`, `constraint_fold_from_source`; saved `KINETIC`, `VIRTUAL_CONSTRAINT_SOURCE`, `THETA_SOLUTION` | Determine which variables the supplied construction actually carries. Do not repeat variation, elimination or the old transcription residual. |
| Completed closed face contributions | c2 `build_case`/`main`; `CLOSED_SLAB_OPERATOR_TERM_ORIGINS` and `CLOSED_SLAB_OPERATOR_PARITY_BLOCKS`, each case's `E_W` entries | Reuse the post-closure, retained-shape, physical-field row returns and already computed half-sum/half-difference combinations. These are pre-`extract` rows, not weakly extracted coupling kernels. Do not recompute them. |
| Closure provenance | c2 `build_face`; each `FOLD_SYMBOL_MAP`'s `REFERENCE_PRESSURE`, `NORMAL_JET`, `PRESSURE`, `DENSITY_BINDING`, `IDENTIFICATIONS` | Establish pressure-versus-reference-trace slots, global normal-derivative orientation, density binding, original source and common closure operator. This is not permission to rerun the kernel/trace map. |

The saved `FOLD_SYMBOL_MAP` values precede the final row restriction. Their
`MULTIGRADE` is the **whole record's aggregate support**, not a per-field
grade declaration. The saved per-face row returns already follow the native
retained-shape projection and physical-field map, before `extract`. Keep those
facts distinct. Preserve the actual amplitude,
independent background/gradient grades, units, profile and coordinate/branch
context; no common scalar or case count authorizes reuse across cases.

## One missing algebraic diagnostic, using completed closures

Use the saved closed face returns rather than re-expanding their pressure
kernels. The following is an implementation route to be checked, not an
assumed proportionality or an expected residual value.

1. Restore the selected original values faithfully with their native symbol
   assumptions and typed metadata. Source/code is read or compiled only as
   needed for faithful decoding; do not import a producer that runs work.
   Persist the exact input references before any new algebra.
2. For each face, extract from the saved centre row and the saved **normalized**
   thickness face-work row the coefficients of that face's pressure and
   global normal-derivative slots. This is new algebra on those complete
   supplied operands, not a repeat of their virtual-work derivation. Save
   the full rows, extraction arguments, coefficients and reconstruction
   residuals before inspecting them. Persist all terms outside those slots.
   A non-slot remainder of the thickness row does not veto its slot-image
   reuse: c2's per-face origin explicitly excludes that remainder. An unjoined
   remainder of the centre row must remain separately reported and prevents
   treating the slot diagnostic as the complete centre load.
3. Determine from those coefficients whether a scalar multiplier **for each face** maps the
   actual normalized pressure/jet face load to the centre load. Derive a
   candidate from the coefficients, then check **both** slot relations and
   the full pressure-dependent rows. Record denominator restrictions without
   dividing by a physical field or silently excluding a parameter branch.
   Test proportionality on the pressure/jet content, not on the entire
   thickness row including its non-slot remainder. If a safe relation is unavailable, stop; do not construct a fresh c2
   closure as a fallback.
4. Reuse the completed c2 face returns only after establishing that this
   face's multiplier can pass through the original closure substitution,
   profile `xreplace`, `retained_shape` and `physical_fields`. There is no
   `extract` operation here. For this bounded route require a scalar independent of
   the wave fields, integration variables, coordinates, profile fields and
   retained small parameters. Check its units against the saved row units.
   Join the actual b `slab_operator/E_W_BALANCE/EXPANDED` pressure-slot terms
   to the normalized face-work row. Pin the b export consumed by c2's actual
   successful manifest and compare its literal `slab_operator` value with
   the b output record; a current filename alone is insufficient. Include
   the saved row-normalization operand. The centre row remains in its original
   action orientation; do not apply that multiplier a second time. If the
   premise fails, record the obstruction and stop.
5. Construct the new centre-load diagnostic from the computed multipliers
   and the **saved** c2 parity combinations, with the source-defined sum and
   difference conventions: with source-defined `S=(Fplus+Fminus)/2` and
   `D=(Fplus-Fminus)/2`, assemble `(rplus+rminus)*S + (rplus-rminus)*D`.
   The two multipliers are independently derived; there is no requirement
   that they be equal and no common-scalar test gates this diagnostic.
   Do not recalculate the old face sum/difference.
   Persist each new operation's full arguments and return, followed by the
   case result, actual support and any zero/nonzero/undecided disposition
   computed from it. Retain formal integrals, outgoing-domain qualifications
   and density/source identities. Do not evaluate new integrals or modes.

The test does not assume a zero result and a nonzero result is not a worker
failure. A failed reuse premise is also a legitimate informative stop. An
undecided simplification must stay undecided; no unbounded simplification
campaign or numerical specialization is a substitute for the generic object.

Two bounded routing controls act at this diagnostic's imported operand boundary:
reverse one pressure-slot contribution, and exchange one face's normal-jet
slot with the other face's slot while retaining its original context. Save
the altered input and actual coefficient/reconstruction response. They test
the new routing and relation checks; they do **not** independently validate
the old acoustic solver or constitute an action-derived physics review.
Do not mutate completed c2 values or rederive them. If a control is degenerate
on the actual input, report that rather than manufacturing a response.

Add one selected final-assembly control: reverse **both** pressure/jet
contributions of one face in the imported centre row, then run the same
coefficient/relation/assembly path. Save its actual multipliers and complete
result even if a routing control stopped earlier. This tests the final parity
assembly as well as the local join. Do not type a predicted control result,
alter a saved parity return, or claim an independent acoustic derivation.
Any zero result inherited through a saved parity return is identified as such
and carries c2's face-sign and other upstream debts. For units, use the actual
nonzero face-row and centre-row dimensions; a stored zero's `[0,0,0]` annotation
does not set the dimension of the physical pairing.

There is no second independent derivation of the closed centre load in this
job. Import/reconstruction equalities and readback are plumbing checks, not
independent physical confirmation. Applicable result/instrument review and
upstream cross-engine debts remain visible; this diagnostic cannot close them.
The [review adjudication](../_measurements/S11c_d_centre_plan_adjudication.md)
retains both original NEEDS REVISION verdicts and records these bounded
corrections. Local correction closure is not a new independent CLEAR verdict.

## Build mechanics and limits

Planned helper: `_measurements/S11c_d_centre_compatibility.py`.
Fresh run: `_scratch/s11c/s11c-d-centre-compatibility-20260926/production`.
Do not populate that run until the focused build-plan gate is adjudicated.
Build → execute the bounded diagnostic → emit/save → report → stop. The
builder does not launch a reviewer or downstream job.

> The script may PRINT computed objects. It may NOT state conclusions.
> PRINT the residual before a guard; do not conceal it behind an assertion.
> Interpretation belongs in the report, not in a typed scientific payload.

The only hand-combined physical expressions are the supplied action/ansatz.
All other physical expressions must be reached by computation from the
saved operands. Controls enter those operands, never a result. Emit every
case irrespective of its value. Do not type a predicted multiplier, centre
equation, parity value or unavailable result as zero.

One guarded algebraic worker only: 900 s whole job, 2 GiB, zero swap, one CPU,
nice 15, 32 tasks, one native-library thread, with the native supervisor
inside `scripts/s11c_guarded_run.py`. No response solve, old producer or
completed scientific operation runs. No automatic retry, overlapping worker,
old-directory restart, unguarded fallback or model polling. Preserve partial
inputs/values and stop at the cap. The silent local hook reports completion
or error. Resource caps are a ceiling, not a runtime prediction.

Acceptance of the **diagnostic record**, not of a centre reduction, requires
actual final guard/supervisor/child exits, strict stderr inspection,
checks/stdout identity, consumed source/input/logical/canonical posthashes and
resource telemetry. Preserve any non-clean run as evidence. A report maps
every result to its inputs/new operations and states the resulting dependency
boundary. Existing real-frequency results, deferred pole work, failure history
and the retained builder-report suffix remain immutable.

## Gate scope

Review this finite diagnostic and its source-defined reuse premises. A
substantive objection must identify a wrong physical input/operation, an
unjustified reuse/restriction, an overclaim, or an inadequacy that prevents
the stated diagnostic. Optional presentation improvements are nonblocking.
Do not reopen the cleared amendment, commission a new centre/radiation solver,
or demand a complete scattering proof from a necessary-condition check.
If the proposed diagnostic cannot answer anything useful within these inputs
and bounds, say so with the specific missing dependency; that is an acceptable
stop, not a reason to keep expanding this plan.
