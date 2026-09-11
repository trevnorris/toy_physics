# S10 coverage and remaining work

Updated 2026-09-11. This is a mathematical and implementation coverage map,
not a claim that the supplied physical model has been derived.

| Requirement | Current coverage | Remaining obligation |
|---|---|---|
| Supplied action and actual first variation | Baseline plus all five controls in Lean | Physical justification remains a premise. |
| Compact-test stationarity and local PDE | All six action families | Boundary/interface and weak-solution extensions are outside this S10 setting. |
| Real cosine phase average | All six actions; exact 0..2pi integral | CAS route-factor/comparator bookkeeping remains separate. |
| Baseline and form-control spectrum | Complete real plane-wave amplitude classification in arbitrary D | This does not prove completeness of general PDE solutions. |
| Anisotropic spectrum and exceptional directions | Complete oblique, parallel, perpendicular, and static classifications in Lean | The focused CAS samples do not certify a general stratum-discovery algorithm. |
| Coefficient and sign controls | Arbitrary real squared-frequency classification; positive/negative sign and N2/N3 counts | An exponential spacetime solution is not yet constructed. |
| Q5 scaling | Nonzero root formulas and ratios, including the anisotropic extra branch | Zero-root ratio remains undefined, rather than assigned a value. |
| Q6 dimensions | PhysLean expression checker and unit-change theorem; actual-action bridges, coefficient solve/inventory, full root ratios, Q7, matrices, mixed-unit minors, complete basis families and N5/N6 residuals; kernel/rank unit invariance and vacuity proved | Connect actual CAS expressions and emissions, including basis normalizations, denominator domains and root substitutions. |
| Q7 ordinary three-dimensional curl | PhysLean Levi-Civita contraction and all six package comparisons, with exact action-stiffness linkage | Align the CAS implementations and production comparator with the explicit construction. |
| Q8 focused CAS repair | Both engines inspect N2 and N3 matrices; ten sampled roots compare successfully | Full production rerun and integration with the existing broad comparator/export chain. |
| Ledger/paper alignment | New evidence recorded in S10 and linked proof reports | Reconcile the paper's historical claims and evidence pointers with the completed formal coverage. |

The focused CAS rerun and its comparison are documented in
[S10_anisotropic_strata_report.md](../../_measurements/S10_anisotropic_strata_report.md).
The latest control results are in [SCALAR_RESULT.md](SCALAR_RESULT.md).
The dimensional-analysis and curl extension is in [Q6_Q7_RESULT.md](Q6_Q7_RESULT.md).
The matrix, minor and complete-basis extension is in [MATRIX_RESULT.md](MATRIX_RESULT.md).

S11c-d currently consumes the frozen S11c-b base and S11c-c1/c2 deltas.
This increment preserves those files and `S10_exports.py`; it makes no changes
to the S11c-d constant-end resolvent or scattering construction. Any proposed
change to those shared physical operands requires a separate dependency review.
