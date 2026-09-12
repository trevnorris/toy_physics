# S10 coverage and remaining work

Updated 2026-09-11. This is a mathematical and implementation coverage map,
not a claim that the supplied physical model has been derived.

## Current scope and Lean completion

[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md) governs this work. The
systematic CAS bridge expansion is stopped; the committed bridge remains
evidence. Outstanding transcript fields are not automatically Lean obligations.

The finite Lean resume plan is:

1. Consolidate the existing six-family proofs into a compact contract identifying
   the action/operator, conventions, domain, exhaustive cases, roots and
   full/transverse kernel counts. Reuse general theorems across dimensions.
2. Check the compact connection to the actual CAS objects and map the relevant
   mutation controls to each load-bearing contract claim. Add work only for an
   identified mathematical or fidelity gap in that contract.
3. Complete the two independent non-author fidelity reviews, resolve findings,
   and record the contract's completion against the build and mutation evidence.

Contract consolidation and review remain open; this policy change does not
declare them complete. Once these items are satisfied, stop the Lean task.
The CAS production, comparator/export and paper tasks below have separate
completion criteria. A discrepancy that invalidates the compact fidelity link
still blocks that link, even when the implementation repair belongs elsewhere.

## Existing coverage and remaining obligations

| Requirement | Current coverage | Remaining obligation |
|---|---|---|
| Supplied action and actual first variation | Baseline plus all five controls in Lean | Physical justification remains a premise. |
| Compact-test stationarity and local PDE | All six action families | Boundary/interface and weak-solution extensions are outside this S10 setting. |
| Real cosine phase average | All six actions; exact 0..2pi integral; actual anisotropic D3 matrices now prove `M_A = -2 M_B` and their action linkage | Retain a compact action/operator normalization link in the contract. Production route comparison remains separate; do not replicate the full expression bridge. |
| Baseline and form-control spectrum | Complete real plane-wave amplitude classification in arbitrary D | This does not prove completeness of general PDE solutions. |
| Anisotropic spectrum and exceptional directions | Complete oblique, parallel, perpendicular, and static classifications in Lean | The focused CAS samples do not certify a general stratum-discovery algorithm. |
| Coefficient and sign controls | Arbitrary real squared-frequency classification; positive/negative sign and N2/N3 counts | An exponential spacetime solution is not yet constructed. That construction is outside the current claim and requires an explicit contract extension. |
| Q5 scaling | Nonzero root formulas and ratios, including the anisotropic extra branch | Zero-root ratio remains undefined, rather than assigned a value. |
| Q6 dimensions | PhysLean expression checker and unit-change theorem; actual-action bridges, coefficient solve/inventory, full root ratios, Q7, matrices, mixed-unit minors, complete basis families and N5/N6 residuals; kernel/rank unit invariance and vacuity proved. Both engines' generic and exceptional anisotropic D3 matrices, determinants, roots, stacks, bases, residuals and emitted minors have checked expression bridges (916 scalar expressions including generic/exceptional counts, root lists, aggregate coincidence differences and spectrum solve operands), explicit domains and complete-kernel basis proofs. Fixed-point reruns have an explicit coordinate/physical scale convention. | Include the necessary dimensional/scaling invariants and parameter map in the compact contract. Remaining per-output metadata and other emitted cases belong to CAS/comparator coverage; no systematic Lean expansion is planned. |
| Q7 ordinary three-dimensional curl | PhysLean Levi-Civita contraction and all six package comparisons, with exact action-stiffness linkage | Align the CAS implementations and production comparator with the explicit construction. |
| Q8 focused CAS repair | Both engines inspect N2 and N3 matrices; ten sampled roots compare successfully. Lean now checks all 83 emitted D3 minors, completeness of the row/column selections, all 12 rank-drop locus predicates and all four targeted points, including both branches of the extra transverse locus. The rerun matrices, roots and all 12 printed basis vectors are connected to complete kernel proofs and physical rescaling. All 112 generic and exceptional rank/nullity/basis-count records are certified against the actual objects; generic bindings keep explicit chart assumptions and N4/N7 use signed subtraction. All six solution/distinct-root list pairs have determinant-completeness proofs; distinct counts, algebraic multiplicities and Wolfram candidate-filter counts are certified. The primary coincidence equations, guarded loci, allowed regions, Boolean/outcome decisions and SymPy witnesses are certified on the positive-coefficient domain, including the full parallel axis. All Q8/aggregate coincidence fields, all 16 reported root signs, eight empty root solver-condition lists, three spectrum solve operands/statuses and six retained/skipped-stratum records are now bound to their mathematical meanings. Lean resolves both SymPy undecided extra-root signs as positive. | Consolidate the proved exhaustive cases and invariants into the CAS coverage contract. Production stratum handling and full comparator/export integration remain separate tasks; further metadata translation into Lean is not required. |
| Ledger/paper alignment | New evidence recorded in S10 and linked proof reports | Reconcile the paper's historical claims and evidence pointers with the completed formal coverage. |

The focused CAS rerun and its comparison are documented in
[S10_anisotropic_strata_report.md](../../_measurements/S10_anisotropic_strata_report.md).
The latest control results are in [SCALAR_RESULT.md](SCALAR_RESULT.md).
The dimensional-analysis and curl extension is in [Q6_Q7_RESULT.md](Q6_Q7_RESULT.md).
The matrix, minor and complete-basis extension is in [MATRIX_RESULT.md](MATRIX_RESULT.md).
The first actual-expression bridge is in [CAS_BRIDGE_RESULT.md](CAS_BRIDGE_RESULT.md),
with current build evidence in [CAS_BRIDGE_VERIFICATION.txt](CAS_BRIDGE_VERIFICATION.txt).
The complete minor and exceptional-locus bridge is in [MINOR_LOCUS_RESULT.md](MINOR_LOCUS_RESULT.md).
The exceptional matrices and complete bases are in [EXCEPTIONAL_RERUN_RESULT.md](EXCEPTIONAL_RERUN_RESULT.md).
The exceptional rank/nullity and signed-count bindings are in [COUNT_RESULT.md](COUNT_RESULT.md).
Root-list completeness and multiplicities are in [ROOT_RESULT.md](ROOT_RESULT.md).
The primary coincidence loci, decisions and witnesses are in [COINCIDENCE_RESULT.md](COINCIDENCE_RESULT.md).
The aggregate/Q8 fields, root signs, solve operands and stratum dispositions are in [METADATA_RESULT.md](METADATA_RESULT.md).

S11c-d currently consumes the frozen S11c-b base and S11c-c1/c2 deltas.
This increment preserves those files and `S10_exports.py`; it makes no changes
to the S11c-d constant-end resolvent or scattering construction. Any proposed
change to those shared physical operands requires a separate dependency review.
