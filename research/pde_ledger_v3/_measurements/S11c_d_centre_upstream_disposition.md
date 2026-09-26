# Centre motion: the remaining upstream decision

Status: **source inspection complete; proposed next scope pending**. This
continues the user-requested inspection after the accepted bounded
[face-load diagnostic](S11c_d_centre_compatibility_report.md). It makes no
new physics claim, changes no governing premise and launches no calculation
or review. The existing Option B deferrals remain.

## What the additional inspection resolves

The [navigation record](S11c_d_centre_upstream_navigation.json) pins the
actual source definitions, selected saved output records and input routes.
It used source ASTs and literal text only: no scientific imports, packet
restoration, algebra, old producer, response solve or new validator.

| Available object | Exact source / saved evidence | What it does and does not supply |
|---|---|---|
| Independent centre direction | a `dof_fields`, both face-source constructors, and the four already indexed `ZETA_C` face velocities | Supplies the independent geometric direction. Material records also contain shared in-plane profile advection; adding the thickness and centre records wholesale would double it. |
| Centre virtual-work test | a `virtual_work_cases`; all four cases' saved `DELTA_W/ZETA_C` and `ZETA_C/ZETA_C` records | Their open pressure/jet work text is byte-identical for the two wave-direction labels. This is a saved-input identity, not a closed centre equation: pressure is still an independent input at this stage. |
| General prescribed-drive acoustic response | c1 `response_operator_case` and `closed_coefficients`; literal `FACE_RESPONSE` and all eight `FACE_RESPONSE_COEFFS` cases | Actual input slots are `V_S` and `MU_THETA`. These are reusable response operands, not a rule selecting the centre velocity. No old face solve needs to be reconstructed for packaging. |
| Restricted parity evidence | c1 `DTN_BY_PARITY`, both anchorings | Both saved `OFF_DIAGONAL_BLOCKS` entries are literal `Integer(0)`. The native construction uses the first-shape bare DtN kernels. This is not a whole closed slab/constraint/current invariant-sector result, a centre equation or a uniqueness proof. |
| Acoustic regularity qualifications | c1 `NONINVERTIBILITY_CONDITION`, four face/anchoring records | These retain a formal operator noninvertibility condition and reserve the profile-conditioned resolvent. An invertible face response cannot simply be inferred from the coefficient table. Nor would it establish invertibility of a separate centre balance. |
| Supplied slab balance construction | b `kinetic_balance_from_energy`, `constraint_fold_from_source`, `face_generalized_force_rows` | The kinetic caller differentiates with respect to in-plane velocities and thickness velocity. The face-row caller uses the thickness wave direction with an independent centre virtual test. These are the inspected constructions; their omission of a centre kinetic row is not authority to add a mass or to declare a complete massless constraint. |
| Existing closed scattering source | c2 `build_face`, `expanded_rows`, `build_case` | Uses `REPRESENTATION='DELTA_W'` and closes `U`, `THETA`, `E_W`. The completed diagnostic establishes that this source produces zero centre load through the saved closure. It does not settle the independent homogeneous sector. |

The earlier uniform S11b construction is not a hidden solution of this gap.
Its geometry specifies the centred thickness displacement; its centre-parity
response is a bulk response to prescribed face motion. The later b/c1
specifications expressly keep the independent face coordinate. The inherited
“no incoming waves from infinity” condition concerns the exterior bulk; by
itself it is not an initial/incident condition for an independent slab-centre
solution.

## Recommended bounded next task

**Complete the upstream centre balance from the existing supplied model,
without adding constitutive terms or prescribing centre zero.** This is an
upstream b/c1/c2 disposition, not an extra mode inserted into d. Approval of
the task would authorize investigating this route, not asserting that it
closes or that its answer is zero.

1. State explicitly how the inherited balance-law/virtual-displacement rule
   applies to the independent centre coordinate. Check whether the supplied
   energy, constraint and external-work data really determine its balance.
   Do not silently interpret a missing term as a supplied zero. If an
   additional kinetic/restoring/holding input is needed, name it and stop.
2. If that model statement is sufficient, specify the centre-dependent face
   contribution from the saved independent direction and generic c1 response
   operands. Retain both physical faces, shifted pressure/global-normal-jet
   conventions, actual density, units and independent grades. Separate the
   shared material-advection term. No repeated old kernel, derivative,
   closure, Taylor, normalization or face-response construction.
3. State the incident/initial and domain conditions that select the response,
   distinguishing them from exterior radiation conditions. Establish only
   the sector restriction actually supported by the complete balance and
   its reverse coupling. A zero forced centre load or a bare-DtN parity zero
   alone does not close this step. If proving selection requires a new
   profile-wide inverse, new mechanics or a solver campaign, stop with that
   concrete requirement.
4. Return one usable upstream relation with its conditions, or one precise
   unresolved model input. Only then adjust the d reconstruction plan.
   Existing five-field results remain preserved on their stated premises;
   neither correctness nor invalidity of their general physical interpretation
   is assumed here.

The task must not become another series of prerequisite diagnostics. Its
first deliverable is the model/balance disposition, not numerical data.
**Rough cost:** small source/specification work to settle the supplied law;
medium implementation effort if a finite missing composition suffices; a
general centre inverse or new mechanics is large/uncertain and outside this
proposal. Any actual missing algebra requires a concrete build/instrument
disposition and the existing single guarded worker. No reviewer is sent a
packet automatically, and no new amendment round is scheduled.

The alternative is to prescribe the centre trajectory as an external input
and label the resulting restricted scattering problem. That is an explicit
scope/model choice, not a harmless implementation shortcut. In particular,
choosing a fixed centre here is not justified by the diagnostic zero and
cannot silently satisfy the current independent-centre requirement.

## Effect on the remaining work

General physical-face reconstruction remains conditional on this upstream
disposition. Bindable representations of the already computed five-field
problem can be prepared, but cannot be relabelled as the complete general
A9/A11/A12 deliverable. The separate radiating-boundary/continuous-spectrum
coverage finding also remains. Completing centre mechanics would not itself
authorize that method or establish bulk escape.

The accepted amendment §2 says: “If no upstream map determines a required
face drive, record an upstream c1/c2 dependency finding and its implications
for any affected end-channel result; do not close the gap by adding or fixing
a mode in d.” This is why a specific upstream disposition is needed before
new physical construction, rather than another row batch or an assumed zero.

The follow-up navigation retained one metadata-only display error: it first
treated the three-field `OFF_DIAGONAL_BLOCKS` tuple as a key/value pair. The
corrected literal navigation handles the label and both saved values. No
scientific operation or file mutation occurred in that failed display.
All selected source/metadata and literal-record posthashes passed. The
navigation report is evidence of locations and saved bytes, not an independent
mathematical review. All old values, failed histories and the protected
builder suffix remain unchanged.
