# S11c-d wider-box quadrature refinement

Production completed cleanly in 10m28s with four single-thread workers. All
40 one-momentum and 30 two-momentum rows were refined on both fields at fixed
cutoffs 3/4, regulator 0.2, position bounds 48/14 and source/profile orders
256/512. The ten three-momentum rows and their 144/24/24 rules remain held.

The final single outer-order change is 2.87e-15 in raw integrals and 2.84e-15
in complete actions. Final paired outer/inner changes are 4.84e-16 / 3.40e-16
in raw integrals. The refined cutoff 3-to-4 complete-action difference remains
1.134674e-7. These sampled refinements are much smaller than that box effect;
independent outer quadrature and wider-box three-momentum checks remain next.

Saved-operand acceptance verifies all 81 current/frozen sources, 28 records,
330 worker artifacts including 274 partials, 40 exact read-only cache records,
6720 zero held-term scalars, every native contraction and literal refinement
array. Production evaluated 4,674,896 new nodes; peak worker RSS was 143,720 KiB.
All workers and the supervisor exited zero with empty stderr; checks/stdout,
pre/post packets and 236,526 metadata paths join. The preflight accepted at
09fb8752 remains instrument evidence only.

The 10,635,882-byte transcript has SHA256
`e6ab9b3c8919412dde2fcc3a9274645575edcf03e9e9e856434e11db3a998e54`.
Its canonical path is `scripts/out/S11c_d_wide_momentum_refinement.out`.
The checkpoint retains source, worker, publication and validation evidence.
Full wider-box convergence, independent-grade/uniform/global exceptional
coverage, infinite tails, Abel limits, scattering and poles remain open.
The engine and approved inputs are unchanged; Q9-only edits have no identified
consumed-input dependency, and pinned exports are not rebased here.

Publication commit `812864e3` stores the canonical transcript through
DataLad/git-annex. Its full SHA256, MD5E key, symlink and Git mode 120000
are verified against the accepted payload.
