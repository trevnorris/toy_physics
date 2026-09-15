# S11c-d two-momentum refinement

The complete run passed in 308 seconds with exit zero and empty stderr, using
one native numerical thread. All 30 two-momentum rows, 30 bound source records
and two nested profile types cover both approved fields and three positions.
The original outer-64/panel-16/source-128/profile-128 baseline replays exactly.
Every held single/three-momentum term is unchanged; actual measure controls respond.

Final raw-integral changes under outer, inner-panel, source and profile refinement
are respectively 8.33e-14, 1.50e-12, 2.62e-16 and 8.55e-15. Adaptive Gauss-Kronrod
outer quadrature agrees within 2.45e-14 in integral arrays and 6.21e-17 in complete
actions, with error estimates at most 1.60e-12. Independence is restricted to the
outer rule on the same inner/source/profile operands. The final action correction
from the old pair baseline is 3.63e-11 / 1.53e-11 for the two fields.

All 22 records, nine exact read-only caches, 20,778 tags, 10,387 fresh keys and
91,047 metadata paths pass. Source and profile cache units remain distinct.
All 52 source hashes, 80 partial and 672 adaptive-point packets are verified;
pre/post-emission packet identities agree. Peak process RSS was 323.4 MiB.
The 5,651,996-byte transcript is published through DataLad/git-annex at
`scripts/out/S11c_d_pair_momentum.out`; SHA256
`9b901f3aac1a46869cf59c5ef0ac5d30c98710f573ee441da18b669fe4fa80f0`.

Next: refine all ten three-momentum raw integrals. Their small weighted action
contribution does not establish raw-integral or independent-grade convergence.
Physical domain tails/interchange, Abel weak limits, matching, scattering and
bound poles remain separate work. Approved inputs and the retained contract stay
unchanged; all accepted operands and failed-run evidence remain in repository scratch.

Publication commit: `3d4a32dd6bc8a6dd6642734d9338df8f7fe132c6`. Full SHA256, annex key, symlink and Git mode 120000 verified.
