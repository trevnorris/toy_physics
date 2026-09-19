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

The complete source transcript is now published and annex-verified: 3,357,686
bytes, SHA256 `e38190f011adac82e4221b931b178d6ceda4f589415399f37eb5422c72055f29`.
Final validation exited zero with empty stderr in 381.61 seconds. All 100
cell reconstructions, 931 independent coefficient derivative scalars and 100
responding coefficient mutations pass. All 43 sources, 127 inputs, 130 artifacts,
686 tags, 341 write keys and 4595 metadata paths validate. The original 170
decoded prefix payloads are identical. All 28 original source packets and 75
inherited proof packets are unchanged; only 25 cell proofs were newly completed.

The earlier guards incorrectly required representation statistics to stay fixed
across serialization/calculation contexts. Live and reconstructed statistics are
preserved; only DAG-node and distinct-derivative counts may vary. All other
physical census fields and exact semantic coefficient proofs remain mandatory.
The last-case snapshot in the successful run matched its earlier restored
snapshot; that does not make these statistics session-invariant. The independent
acceptance checker normalized the saved mutation operand with the original
`expand` operation; its initial tree comparison and logs are preserved.

The first recovery completed 75 cells, 698 independent coefficient derivatives
and 75 responding coefficient mutations, plus those cases' source/limit checks.
The final validator reuses that work, validates the last 25 cells, and completes
the unchanged output and metadata replay. The narrow scheduling/statistics
adapter has a reverse whole-function AST join and rejects all 12 tested physical
census mutations. No source construction, integration or solve is repeated.

Use the exact matches to reuse numerical operands for the remaining responses. Their continuum,
current, control, scoped pole and final export work remains. The selected
32-point pole diagnostic is closed without a resolved candidate; no additional
angular doubling or broad quadrature campaign is queued.
