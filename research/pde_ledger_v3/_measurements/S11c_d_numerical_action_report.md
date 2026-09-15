# S11c-d numerical action checkpoint

The user approved the proposed values for all 30 inherited gradient-energy
coefficients. Commit `93692321` records the separate development input and its
executable projection check: every existing parameter, profile and reference
unit is unchanged. The exact dependency census covers both strong and weak
pencils at LEFT, RIGHT and REFERENCE; none contains the added coefficients.
The previous end and normalization records retain their original input hashes.

The existing engine now contains NumericalReducedAction for binding the saved
operator and BoundedActionQuadrature for evaluating its nested integrals. The
adapter keeps the four local derivative orders and all saved nonlocal operands.
It substitutes two Gaussian test ansatzes directly into the full reduced rows
and separately applies the extracted local coefficients and nonlocal operators,
covering all five input fields, all five equations and three positions.

The numerical checker preserves every contribution, both action routes and
their residuals on two finite quadrature grids. The routes contract weights
using separate sum and dot implementations. Native limit order, all momentum
roles, measures, profile subtraction and finite Abel regulator are explicit.
It also perturbs an actual local coefficient and an integration measure to
measure sensitivity. Operands are saved before emission and after each grid.
Full source, packet, dimension/grade and emission replay guards precede acceptance.

The quadrature smoke check against an independently integrated bounded nested
Gaussian agreed to 8.18e-14. The native numerical run is accepted on its finite numerical domain.
Finite shared-grid agreement tests implementation; spatial/momentum tails,
quadrature convergence and the Abel weak limit must still be resolved before
two-ended boundary matching. No S-matrix or profile-frequency pole is produced
by this checkpoint. Existing symbolic operators and the builder report's
retained solver/export contract remain unchanged.

Implementation is committed at `2663a526`. The native constructor and validator
completed with exit zero and empty stderr in 649.93 seconds (167,712 KiB peak
RSS). All 300 field/equation action evaluations agree within 2.49e-16; all
1,920 nonlocal contributions are retained. Each route evaluated 86 distinct
integrals, including the six nested profile transforms. Coefficient mutations
changed the selected contribution by 0.0037616, and measure mutations by
0.0002133–0.0002511. Full replay covers 106 tags, 3,611 metadata paths and
51 fresh write-keys; the packet hash is unchanged across emission.

The two grids change the actions by up to 0.03216 and 0.04670 for the two test
fields. These differences require quadrature refinement; no physical convergence
is claimed. Source/input/packet hashes were checked again before publication.
The 287,329-byte transcript is published as scripts/out/S11c_d_numerical_action.out
through DataLad/git-annex, with SHA-256
3855a4300c4c9a5116b70697e32678e1258a270dbebc70357d7195c2599e6159.
