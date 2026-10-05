# S10 D3 aggregate records, root signs and solver metadata

This extension adds 48 metadata records and nine arithmetic records containing
54 scalar expressions from the focused anisotropic D3 transcripts. Six Wolfram
aggregate records receive both checks, giving **51 newly selected transcript
tags**. The combined arithmetic bridge now covers **289 records and 916 scalar
expressions**.

| Selected metadata | Records | Checked meaning |
|---|---:|---|
| SymPy Q8 copies of coincidence fields | 9 | Each copied locus, allowed operand and Boolean decision agrees with its primary record and proved geometry |
| Wolfram coincidence associations | 6 | Every pair, equation, locus, allowed region, outcome and Boolean decision is checked; the pair list is complete and ordered |
| Root sign observations | 16 | The reported observation is preserved and the actual imported root has an independent sign proof |
| Wolfram root solver-condition lists | 8 | The printed residual-condition list is empty; the imported value is a determinant root on the explicit coefficient domain |
| SymPy spectrum solve statuses | 3 | The returned candidate list is nonempty and exactly covers the zeros of the printed solve operand |
| SymPy retained/skipped stratum records | 6 | Both retained targets are admissible; the skipped origin conflicts with positive wavevector norm |

The arithmetic records comprise 48 scalar entries in the Wolfram associations
(two indices and one root difference for each of 16 pair records), plus the
three SymPy solve-operand/solve-variable pairs. Their expressions have checked
values and units. Original payloads, line numbers, hashes, trees, schema fields,
and theorem names are retained in the
manifest (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_manifest.json`).

## Complete aggregate and duplicate-field bindings

RecordBindings.lean (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/S10Audit/CAS/RecordBindings.lean`) translates each repeated
field from its own payload. It proves an identity with the primary predicate,
then transfers the primary geometric theorem. Each association's printed
root difference vanishes exactly on its printed locus. Its Boolean and outcome
decisions agree with existence in its own allowed region. Pair-index lists are
proved equal to `[(1,2),(1,3),(2,3)]`, or `[(1,2)]` for the parallel rerun's two
distinct roots.

The association parser requires exactly the six declared keys, each once.
It preserves the list order and all guarded branches; duplicate or missing
keys are rejected. Equations must compare to literal zero, and pair indices
must be integers in 1..3. Wrong in-range indices and wrong root-difference
coefficients are rejected by Lean's mathematical proofs. These checks cover
the generic aggregate, all three Q8 aggregate copies, and both targeted
parameter-locus aggregates.

The coefficient domain and the treatment of excluded square-root branches are
unchanged from [COINCIDENCE_RESULT.md](COINCIDENCE_RESULT.md): D=3,
positive `rho`, `mu`, `sigma`, and `lambdaScale`, with `sigma != 1`.
The coincidence geometry includes arbitrary real wavevectors and the full
exceptional parallel axis.

## Reported signs and mathematical signs

RecordSupport.lean (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/S10Audit/CAS/RecordSupport.lean`) distinguishes zero,
positive, negative, and undecided sign observations. An undecided observation
asserts no sign. It also proves the reference static root is zero and both
nonzero branches are positive when `rho > 0`, `mu > 0`, `sigma > 0`, and
`k != 0`.

StatusBindings.lean (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/S10Audit/CAS/StatusBindings.lean`) applies that theorem to
each **actual imported root expression**, including each engine's own targeted
point. Every sign record has separate theorems for its reported observation,
its computed mathematical sign, and consistency of the reported information.
The computed-sign theorem is required even for an undecided observation.
Generic sign proofs require nonzero wavevector, without a generic basis-chart
restriction.

The transcripts contain six zero observations, eight positive observations,
and two undecided observations. SymPy leaves the extra root undecided in the
generic and perpendicular cases; Wolfram reports both positive. Lean proves
both positive on the declared domain, giving six zero and ten positive
computed signs. The original SymPy observations remain unchanged. This resolves
the mathematical uncertainty in those records without claiming that SymPy
itself decided the signs.

## Solver operands and stratum dispositions

Each SymPy spectrum operand is parsed together with its exact solve variable,
`omegaSquared`. A same-unit composite expression is rejected as a solve
variable. The operand has a checked equality to the reference route-B
determinant, and its complete zero set equals the actual returned candidate
list. The `roots_returned` flag is linked to that list's nonemptiness. These
proofs describe the selected returned data and do not certify the solver
algorithm on arbitrary inputs.

All eight Wolfram residual-condition lists are empty. Each associated imported
root is separately proved to annihilate the reference route-B determinant.
The empty residual lists do not remove the coefficient assumptions used by
these proofs.

The printed skipped branch is proved to be exactly `k = 0`; every point on it
has zero wavevector norm. The primary and aggregate `False` decisions agree
with the impossibility of a positive-norm point on that branch. The printed
reason is therefore justified. The two `not_skipped_allowed_branch` records
are tied to the actual parallel and perpendicular targets: each is nonzero
and belongs to the already certified extra-transverse rank-drop locus.

## Verification

S10_lean_cas_records.py (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/S10_lean_cas_records.py`) extends the
non-evaluating bridge. Independent parser comparisons cover all 48 complete
metadata payloads and all 54 new arithmetic expressions. The comparisons use
the established SymPy parser and the independent Mathematica tokenizer and
full-form parser, retaining the earlier explicit handling of membership atoms
and empty Wolfram lists.

New malformed-input tests reject duplicate/missing association keys, invalid
pair indices, nonzero equation right sides, unsupported outcomes/signs, missing
pairs, and a composite solve variable. New Lean mutations change a Q8 decision,
a pair index, a root-difference coefficient, the extra root's sign, the nonzero
wavevector premise, and the skipped-origin disposition. Positive controls check
a Q8 copy, the preserved undecided observation alongside its positive-root
proof, and complete roots of a printed solve operand.

The full build passes **1962 selected axiom audits**, including **1712 CAS
audits**, across **115 canonical Lean files** and 3834 build jobs. This extension
adds 334 selected audits. Thirty declarations require no axioms; the others use
only `propext`, `Classical.choice` and `Quot.sound`. Warnings are errors and there
are no proof admissions.

The complete regression suite passes 916 arithmetic comparisons, 56 primary
coincidence payload comparisons, 48 metadata payload comparisons, 324 locus
sign-pattern comparisons, four point comparisons and three filter-list checks.
All 80 malformed inputs and 36 mathematical mutations are rejected; all 14
positive controls pass. Canonical source, input and frozen export hashes remain
unchanged by the tests.

Current full-build and regression evidence is in
CAS_BRIDGE_VERIFICATION.txt (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/CAS_BRIDGE_VERIFICATION.txt`) and the
check record (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_checks.json`).
The lexer/parser and transcript translation remain tested software; the
translated mathematical claims are checked by Lean's kernel.

## Remaining S10 integration

The primary, aggregate and Q8 coincidence records, root signs, root solver
condition lists, spectrum solve operands/statuses, and skipped-stratum
records selected here are now connected to their mathematical meanings.
The former plan to translate the remaining D3 metadata and D4/other-package
outputs into Lean is superseded by
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md). Resume with the compact
coverage contract, action/operator fidelity link, mutation mapping and fidelity
review. Remaining CAS metadata coverage, production Q7, the broad comparator/export
chain and paper alignment have separate obligations in [COVERAGE.md](COVERAGE.md).

No CAS engine, input transcript, frozen export or S11 operand is changed by
this extension. The combined scope and resume plan are recorded in
CAS_CHECKPOINT.md (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/CAS_CHECKPOINT.md`).
