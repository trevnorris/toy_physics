# S11c-d quadrature and limit checks

The complete bounded action check is published at `009667f7`. Source and
assembled actions agree to 2.49e-16, but changing the two quadrature grids
changes the action by up to 0.04670 in the declared unit frame. Continue step 5
of the numerical action plan using the saved operands; do not rebuild the
reduced source or the test-field substitutions.

1. Inventory every source denominator that remains coordinate dependent after
   numerical binding, including the three bulk-momentum sites and all three
   Abel transfer pairs. Retain original and bound operands, real/imaginary
   parts and squared moduli. Record positive-regulator and real-momentum
   assumptions; sampled distances alone do not establish global coverage.
2. Transform the already computed Abel even kernel to the physical normal
   momentum coordinate using SymPy's integral change of variables. Derive its
   width and bounded antiderivative from that operand. Compare ordinary Gauss
   grids with the exact bounded mass at both centered and off-node peaks.
   Compare a rule split around the computed width with the same exact mass.
   Include all transfer-pair records; a generic momentum sample misses their
   diagonal concentration. This checks resolution at positive regulator, not
   a pointwise substitute for the distributional weak limit.
3. For every actual nested profile transform, derive its transfer coordinate
   from its phase. Bind only the approved profile instance, retain its original
   leg assignment, and evaluate zero and nonzero transfers. Compare split
   Gauss integration with independent adaptive integration at two spatial
   cutoffs. Record the absolute profile-tail integrals and integration error
   estimates. The zero-jet subtraction remains intact on both half-lines.
4. Use these results to select the next full-action refinement. Bound memory
   before increasing tensor sizes; split or change integration coordinates
   around computed concentration loci when needed. Verify the change of
   variables, Jacobian and source action joins. Reuse all saved direct and
   assembled test operands. Refine spatial/momentum domains and regulator only
   after resolving each fixed-regulator quadrature.
5. Establish the needed action-level tail and Abel weak limits before boundary
   matching. Elementary kernel/profile tests supply numerical operands and
   resolution evidence; they do not establish the full physical operator limit
   or an S-matrix. Retain exceptional domains and failures explicitly.

Preserve the full symbolic input and the builder report contract suffix. Keep
durable runs in repository _scratch/s11c/, publish validated .out files using
DataLad/git-annex, and commit each step. Use one CAS job and a silent local
completion/error watcher, with no recurring checks or model polling.
