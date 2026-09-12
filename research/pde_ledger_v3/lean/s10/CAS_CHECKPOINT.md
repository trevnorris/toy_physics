# S10 CAS bridge checkpoint — 2026-09-11

**Scope update after this checkpoint:** follow
[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md) and
[COVERAGE.md](COVERAGE.md). The bridge expansion plan has been superseded;
the recorded proof and verification evidence below is preserved.

This checkpoint records the S10 work after the
[original S9/S10 checkpoint](../CHECKPOINT.md), commit `56595cf7`.
It adds the bridge from the two focused anisotropic D3 CAS transcripts to
Lean expressions and mathematical proofs, together with its generators,
regression instruments, manifests, verification record and ledger updates.

## Included coverage

- Actual generic and exceptional matrices, determinants, roots, stacks, bases
  and residuals have checked values, dimensions and explicit denominator
  assumptions. The route normalization is `M_A = -2 M_B`.
- All 83 emitted D3 minors, their row/column selections, 12 exceptional-locus
  predicates and four targeted points are checked. The full parallel and
  perpendicular loci are covered, including both vectors at parallel double
  roots. Targeted reruns have an explicit restoration of physical units.
- All 112 generic and exceptional rank/nullity/basis-count records are bound
  to the actual matrices and complete bases. Generic chart restrictions remain
  explicit; count residuals retain signed subtraction.
- Solution and distinct-root lists have determinant-completeness proofs,
  checked algebraic multiplicities and candidate-filter counts.
- Primary, aggregate and Q8 coincidence fields have checked root differences,
  guarded loci, allowed regions, decisions and witnesses.
- All 16 reported root signs, eight empty root-condition lists, three spectrum
  solve operands/statuses and six retained/skipped-stratum records are bound
  to their meanings. Lean proves positive the two extra-root signs left
  undecided in the SymPy transcript, on the stated coefficient domain.

The arithmetic bridge contains 289 tagged records and 916 scalar expressions.
It also checks 56 primary coincidence logical/container payloads and 48 further
metadata payloads. The [manifest](../../_measurements/S10_lean_cas_bridge_manifest.json)
retains exact payloads, source hashes, expression trees, domains and proof names.
Translation remains tested software; Lean's kernel checks the translated
mathematical claims.

## Verification at the checkpoint

The five Lean libraries build successfully with warnings treated as errors:
**1962 selected axiom audits**, including **1712 CAS audits**, across **115
canonical Lean source files** and 3834 build jobs. Thirty selected declarations
require no axioms; the others use only `propext`, `Classical.choice` and
`Quot.sound`. There are no proof admissions or custom axioms.

The [regression record](../../_measurements/S10_lean_cas_bridge_checks.json)
reports 916 independent arithmetic comparisons, 56 primary coincidence and
48 metadata payload comparisons, 324 coordinate sign-pattern comparisons,
four targeted point comparisons and three filter-list checks. All 80 malformed
inputs and 36 mathematical mutations are rejected; all 14 positive controls
pass. Mutation checks preserve canonical sources and inputs.

All 135 artifact hashes and the full build-log hash in
[CAS_BRIDGE_VERIFICATION.txt](CAS_BRIDGE_VERIFICATION.txt) were checked against
the working tree before checkpointing. Reproduction commands are in
[CAS_BRIDGE_RESULT.md](CAS_BRIDGE_RESULT.md).

## Resume point and remaining scope

S9's original formalization pilot remains complete within its declared scope.
S10's mathematical core covers all six supplied action families, but S10 is
not complete end to end. The user subsequently stopped systematic CAS bridge
expansion. The former plan to translate the remaining D3 metadata and then all
D4/other-package emissions is superseded by the scope policy.

The subsequent compact mathematical/coverage contract is now complete: its
action/operator fidelity link and mutation coverage were checked, and fresh
Claude and Grok statement-fidelity reviews both returned CLEAR. See
[FIDELITY_REVIEW.md](FIDELITY_REVIEW.md) for evidence and dispositions and
[COVERAGE.md](COVERAGE.md) for the closed C1–C4 criteria. Lean work stops at that
scope; this historical bridge checkpoint does not create another expansion task.

Production Q7 and stratum handling, the broad comparator/export pipeline and
ledger/paper reconciliation retain their separate implementation obligations.

An explicit exponentially growing spacetime solution for the negative-stiffness
control is also not yet constructed; adding that claim requires an explicit
contract extension. The supplied action and physical dimension remain premises.
See [COVERAGE.md](COVERAGE.md) for the full boundary.

This checkpoint changes no CAS engine, input transcript, frozen export or S11
file. The original full-sweep artifacts and `S10_exports.py` are preserved.
Concurrent S11 work is excluded. Build caches, temporary mutation sources and
transient logs remain ignored; durable verification evidence is included here.
