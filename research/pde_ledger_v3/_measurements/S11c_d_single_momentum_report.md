# S11c-d one-momentum refinement

The complete targeted run passed in 114 seconds with exit zero and empty stderr,
using one native numerical thread. All 40 one-momentum rows and their 40 bound
source records cover both approved fields and three positions. The accepted
outer-64/source-128 rows and complete actions are reproduced exactly. Every
held two/three-momentum term remains unchanged, and actual weight controls respond.

Outer orders 64/96/144/216 reduce the final integral change to 6.64e-15. Source
orders 128/192/256 change integrals by at most 3.81e-15. Independent adaptive
Gauss-Kronrod outer integration agrees within 1.85e-15 in the integral arrays
and 1.90e-15 in complete actions; its largest estimated error is 8.89e-12. This
checks the outer rule on the same source operands, not independent source physics.
The refined action differs from the old finest grid by 1.77e-7 / 4.78e-8.

All 14 numerical records, six exact read-only rule caches, 13,230 tags, 6,613
fresh keys and 64,385 metadata paths pass. Every source/record hash is verified;
all accepted and new numerical packet hashes are unchanged by emission. Peak
process RSS is 261.3 MiB. The 3,838,693-byte transcript is published and
annex-verified at 9fad1d11: `scripts/out/S11c_d_single_momentum.out`; SHA256
`6176e624f8eb852bb7de3bd3a2a03a99592db2dcb3980f4d41c62544b0a8e2cf`.

The one-momentum discrepancy is resolved for these finite numerical tests.
Next: refine the two-momentum contribution, retaining this result and the saved
three-momentum terms explicitly. Held-layout changes, independent parameter
grades, physical domain tails/interchange and Abel weak limits remain distinct
requirements. No matching, scattering or pole result follows from this stage.
All approved inputs, saved operands and the solver/export contract are preserved.
