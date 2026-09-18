# D5 bulk fidelity record

Bounded D5B.1–D5B.4 complete: local verification PASS and both independent
fidelity reviews CLEAR.
Scope is D5B.1–D5B.4 in `D5_BULK_COVERAGE.md`.

## Intended object and theorem boundary

`Coeff = Fin 3 -> R` represents `(a,b,c)` independently of the five field
components. `Point 5 = Fin 6 -> R` has time row zero; spatial row `i.succ`
is native `x_(i+1)`. `spatialGradient J i j = J i.succ j`, so the derivative
index is the row. The actual Lagrangian is
`-1/2 [a (div J)^2 + b sum_ij J_i,j J_j,i + c sum_ij J_i,j^2]`.
`density_identity` connects it to the reviewed D5 invariant form and
`all_invariant_densities` uses the existing unique full SO(5) classification.
The O(5) family is the same by that completed classification.

Momentum is the ordinary derivative of this Lagrangian along a basis jet,
not a supplied tensor. Spatial momentum is
`p_ij=-(a delta_ij div u+b partial_j u_i+c partial_i u_j)` and its time row
vanishes. The local EL is defined as `-sum_j partial_j p_ji`, with component
and derivative indices fixed by the typed function arguments. Coefficients
are constant parameters, not varied fields.

The action is the actual integral of the density change under a smooth
compactly supported variation of a smooth background. Existing S10
integrability, differentiation and integration-by-parts lemmas give its
first derivative and equivalence of stationarity with the local EL equation.
No infinite background action or unjustified differentiation under an
improper integral is assumed. The already dimension-general coordinate rules
in `S11D4Odd.Calculus` supply mixed-derivative symmetry and product rules.

The proved local equation is `(a+b)grad(div u)+c Delta u`. Universal equality
of these operators is equivalent to equal `(c,a+b)`; smooth transverse and
longitudinal plane waves prove necessity while the full local formula proves
sufficiency. Surjectivity of this linear response map gives image dimension
two; its kernel consists of `(t,-t,0)` and has dimension one. Non-nullness
also implies existence of a smooth background and compact test field with
nonzero actual first variation, by the proved stationary/EL equivalence.

The current `J_i=sum_j(u_i partial_j u_j-u_j partial_j u_i)` has divergence
`(div u)^2-tr(G^2)`. A null density therefore equals `-t div J/2`.
Its nonzero affine-field witness is retained: density `-1`, current
divergence `2`. Momentum for the normalized trace-square action on that
jet is `-2`. Zero bulk variation does not remove boundary effects or imply
pointwise zero density. The compact homogeneous map is
`rho=0, mu=c, B=a+b+c` at the bulk/modal operator level, without a pointwise
action-density equality or a new spectrum/stability claim.

## Compact native identification

`_measurements/S11_lean_d5_bulk_source_check.py` selects the original
`compute_q9`, `q9_v5` and coordinate/EL helpers from the existing SymPy source
AST. It does not import the production module or execute its driver. It
compares all three actual native invariant basis vectors with the trace
basis through an invertible exact change of basis, actual `L=-Q/2` momenta,
all V5 responses, modal signs and the homogeneous parameter map. The native
helper uses `+div(momentum)`: for the same Lagrangian its result is the
negative of Lean EL; native V5 applied to Q is twice Lean EL for `-Q/2`.
The report establishes this on all three basis elements.

The instrument also tests original native EL sign/factor mutations, a wrong
transpose in the b momentum, an incorrect homogeneous map, null/non-null
and nonzero normalization witnesses, and responses along the fifth spatial
direction. Wolfram anchors are source inspection only. This is a tested
translation outside Lean's kernel; no CAS software correctness theorem is
claimed.

## Recorded verification and preservation

Five modules plus the audit root expose 54 selected theorem audits. Fourteen
paired false claims are matched with true controls; four additional positives
retain unrestricted null coefficients, zero wavevector, negative coefficients
and existence of nonzero variation. Eighteen positive executions include a
repeated momentum-normalization source for its separate sign/factor pairs.
Each paired rejection has exactly one `contract_control` error and one
unsolved `False`, without warnings, syntax/import/resource failures.

Verification uses isolated output objects and one Lean worker, `-j1 -M4096`,
strict warnings and a 600-second process limit. Timeout cleanup is the tested
portable-runner whole-process-group function. The fresh portable D5 replay
provides the 42 classification dependency objects actually imported here;
its entire 45-object record, inputs, logs, packages, direct source/object
hashes and positive/mutation outcomes were validated first. This is recorded
dependency reuse, not a new execution of those canonical builds. Twelve
other unchanged imports and six new module/root objects are compiled in the
isolated directory. Further reuse requires the exact command, full transitive
source/pin/generator and input/output object hashes. Controls always run fresh.

The preservation manifest pins 499 historical source/evidence files, all 232
existing shared ledger objects, two native source inputs, and the original
installation-validation record. Those records remain historical. Current
lakefile/portable-runner edits are declared D5 integration changes, not silent
updates to prior reviewed packets. The full transitive Mathlib cache remains
the pinned baseline; direct imports and clean package revisions are checked.

## Observed results and provenance

Recorded run4 passes 92 records: 60 canonical object records, 54 standard-axiom
audits, fourteen paired mathematical rejections and eighteen positive
executions (seventeen distinct positive statements). Every false statement
has one `contract_control` diagnostic and one unsolved `False`; every positive
compiles without warnings. No canonical source-replacement mutants are claimed.

Of the sixty objects, 42 unchanged D5 classification dependencies were copied
from the fresh completed portable D5 replay after validation of its full
45-object record. Twelve other unchanged imports plus Action/Variation were
built in recorded bulk run1; Bulk in run2; Census/Controls/root in run3. Run4
reuses those objects under complete guards and executes all 32 controls fresh.
The 54 axiom lists come from the recorded root build in run3, whose source,
command, input/output objects and audit log were revalidated; they were not
printed by a fresh root execution in run4. Only standard Lean axioms occur.

The native check passes twelve exact identities, two mutations of the original
native EL helper (sign and factor), and nine target checks: two wrong computed
formulas, one generic polynomial factor sanity check and six positive witnesses.
The factor sanity check does not independently mutate the computed current;
actual current normalization is supported by the native divergence identity
and Lean current/density witnesses and paired control. Native-to-trace rows are
`(0,1,0),(1/2,-1/2,0),(0,-1,1)`, determinant `-1/2`; target coefficients
`(a,b,c)` have native weights `(a+b+c,2a,c)`. Native run1 evidence is retained
with matching live source/instrument/report hashes, without an unnecessary rerun.

Author validation checks all live sources, sixty objects, archived build
lineage, logs, native evidence, fifteen clean package pins and thirteen direct
Mathlib source/object pairs. It also checks the 499 historical source/evidence
files, 232 shared objects, two native inputs and installation records.
`INSTALL_D5_VALIDATION.json` separately records the full portable D5 replay;
its predecessor `INSTALL_VALIDATION.json` remains historical and unchanged.
No S11c source or pinned export was changed or executed.

## Bounded proof-code repairs

Run1's EL proof hit the default heartbeat limit during repeated coordinate
expansion. Replacing only that proof body with an arbitrary-index momentum
derivative and finite-sum argument passed strict focused and recorded builds.
Run2 then exposed one unused simplification argument in Census; removing it
preserved every definition and statement. Run3 built the remaining modules,
passed the audit and all fourteen mutants, but an extra negative-coefficient
positive stopped at an unreduced third vector entry. Run4 adds
`Matrix.cons_val_two` only to that and the remaining nonzero-variation positive
proof. Both focused cases passed; all canonical sources and objects stayed
unchanged. Full source/object/log/preservation validation preceded each repair.
No resource, warning, positive-control or other instrument failure was counted
as a mathematical rejection. Limits and acceptance/reuse rules were unchanged;
focused checks emitted no object and supplied no recorded build evidence.

## Independent reviews and stopping boundary

The user authorized Claude and Grok review of the current D5B packet and a
commit after both sign off. The verified source/evidence packet was frozen before
transfer. Both independent substantive fidelity reviews returned CLEAR with no required
correction; local compilation and this author's validation are not review legs.
See `D5_BULK_FIDELITY_REVIEW.md` for the fixed revision, reviewer identities,
optional-note dispositions and source-reading limits. The nonzero-first-variation
existence theorem is nonconstructive; no explicit test field is claimed.
All exclusions in `D5_BULK_COVERAGE.md` remain binding. Stop after D5B.1–D5B.4
and their reviews; no broader physics or CAS production completion is claimed.
