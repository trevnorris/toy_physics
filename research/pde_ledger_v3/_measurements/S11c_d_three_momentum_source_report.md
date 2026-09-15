# S11c-d three-momentum source/profile refinement

The initial three-momentum result is published at 22e78b55 and annex-verified
at b87885b5. Its final raw outer/innermost/middle changes are below 2.68e-12 /
8.22e-15 / 5.61e-14. This next stage keeps the accepted momentum grid and all
single/pair contributions, doubles source order to 256, then profile order to
256. It covers all ten rows, both approved fields and all four nested profiles.

The engine is unchanged. The new checker verifies the accepted packet chain,
reconstructs the retained mixed baseline, and emits the actual distinct momentum
orders for every held layout. It saves all complete and partial results and
preserves raw integral changes alongside weighted term and full-action changes.
The accepted baseline is reused explicitly. No source/profile full refinement
result is accepted yet; independent quadrature, parameter-grade coverage,
physical tails/interchange and Abel weak limits remain separate work.

Focused checks passed in 21.5 seconds with exit zero and empty stderr. Both
16,384-node baseline prefixes match exactly in every value, mutation, mass,
node and weight. The complete retained action also matches exactly. Cached
source/profile rules are exact and read-only; actual measure and order controls
respond. Higher-order prefixes are evaluated without claiming full integrals.
The retained-baseline emission replays 2,046 tags, 1,021 keys and 17,637 metadata
paths, with explicit momentum-limit, box-mass and source/profile cache units.

The next four new full grids evaluate 138,624,384 momentum nodes. Measured
prefix rates suggest 2.1-3.4 hours before checkpoint/emission overhead; allow
roughly 3-4 hours. This is a runtime estimate, not convergence evidence.
