# Remaining case source preparation

All three missing native source/assembly cases were saved by 173.24 seconds.
The accepted LAB_HELD/RHO4_CONSTANT case is copied unchanged. The saved new
case counts are:

| Case | Integral operands | Nonlocal cell terms | Exact baseline integral matches |
| --- | ---: | ---: | ---: |
| LAB_HELD / RHOBR_CONSTANT | 70 | 140 | 61 |
| MATERIAL_ADVECTED / RHO4_CONSTANT | 80 | 185 | 61 |
| MATERIAL_ADVECTED / RHOBR_CONSTANT | 70 | 162 | 46 |

These are source counts, not new scattering results. Every coefficient and
ordered integral remains available; an identical integral does not imply an
identical complete operator or response.

Final acceptance is pending output validation. Two guards incorrectly required
expression representation counts to stay fixed across serialization/calculation
contexts. Both live and reconstructed statistics remain preserved. Only DAG-node
and distinct-derivative counts differ; physical census fields remain exact.
Acceptance uses the complete source/coefficient identities and exact semantic
proofs, with no nonzero residual waived.

The first recovery completed 75 cells, 698 independent coefficient derivatives
and 75 responding coefficient mutations, plus those cases' source/limit checks.
The final validator reuses that work, validates the last 25 cells, and completes
the unchanged output and metadata replay. The narrow scheduling/statistics
adapter has a reverse whole-function AST join and rejects all 12 tested physical
census mutations. No source construction, integration or solve is repeated.

After full validation, publish the source transcript and use the exact matches
to reuse numerical operands for the remaining responses. Their continuum,
current, control, scoped pole and final export work remains. The selected
32-point pole diagnostic is closed without a resolved candidate; no additional
angular doubling or broad quadrature campaign is queued.
