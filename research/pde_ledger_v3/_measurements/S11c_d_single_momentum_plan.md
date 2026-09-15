# S11c-d focused one-momentum refinement

The complete finite momentum run is published at 20fe9381 and annex-verified
at 621d81db. Its final full-action changes are 3.64e-5 / 1.05e-5. Decomposition
of the saved native term arrays shows that the one-momentum layout supplies
3.64e-5 / 1.05e-5, the two-momentum layout at most 7.72e-9, and the
three-momentum layout at most 2.23e-12. The decomposition reproduces each full
action change to 2.90e-16. This is an error-location diagnostic, not a proof of
convergence of the held layouts or of individual multigrade coefficients.

1. Verify every accepted momentum source, artifact and saved row/group/partial
   hash. Consume its bound operand packet and finest complete grids directly.
   Restore only the saved reduction/metadata context. Add a SingleMomentum
   helper while proving the accepted engine AST is otherwise identical.
2. Retain the exact 40 one-momentum rows and their 20 distinct source integrals
   for each approved test. Preserve the other 40 native rows at their accepted
   finest quadrature. This intentionally varies one numerical rule at a time;
   every resulting complete action must identify the held and changed rules.
3. Reproduce the accepted outer-order-64/source-order-128 layout first. Evaluate
   outer orders 96, 144 and 216 at source order 128; then source orders 192 and
   256 at outer order 216. Keep every finite bound and the regulator unchanged.
   Reconstruct every native term and full five-field action from the actual
   coefficients, joining all held terms to their accepted values. Save each
   group and complete action before subsequent checks.
4. Independently integrate the one-momentum vector with adaptive Gauss-Kronrod
   quadrature at source order 256. Retain the same literal factor operands;
   independence applies to the outer quadrature, not the integrand or source
   quadrature. Record interval contributions, estimated errors and differences
   for all 40 rows and all three positions. Use absolute/relative tolerance
   1e-10 in the declared numerical unit frame; retain nonconvergence explicitly.
5. Cache only the unchanged one-dimensional source Gauss rules. Compare every
   cached node/weight array exactly with the original rule; make cached arrays
   read-only. Set native numerical thread pools to one for this focused run.
   Reproduce the old layout with the new evaluator before accepting differences.
   Retain all workspace counts separately from process RSS.
6. Emit original/fresh/held operands, all component and term differences, actual
   momentum-measure controls, full units/grades and fresh write keys. Keep
   structural manifests lossless; require full replay and pre/post source and
   packet hashes. Publish accepted evidence through DataLad/git-annex and commit.

The held two/three-momentum terms still have their measured finite-grid changes.
No mixed-grid action establishes uniform momentum convergence, infinite tails,
interchange, an Abel weak limit, boundary matching, scattering or poles. Inspect
the observed differences to select the next numerical step. Keep all approved
inputs and the builder report contract. Use repository scratch, one supervised
job and a silent local completion/error watcher. Never rebuild the accepted
four-hour momentum grids or upstream physics for emission/schema plumbing.
