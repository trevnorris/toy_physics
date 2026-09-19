# One finite-frequency angular refinement

The 16-point contour is accepted at395e4ea9. All matrices are regular and the
8/16 sampled winding is zero, but maximum phase changes remain2.70 radians.
The first inverse moment shrinks from8.74e-2 to7.02e-6; response moments at16
points are below9.24e-14. This warrants one angular comparison before interpreting
the numerical search. It does not warrant a broad frequency or quadrature grid.

Keep the same circle centred at1-0.01i with radius0.02, physical parameters,
finite discretization, positive regulator and approximate modal boundaries.
Reuse all16 accepted matrices and inverses by immutable artifact addresses;
compute only the16 midpoints. Four single-thread workers with2GiB ceilings own
four new points each. Initial wall budget900seconds; the preceding equivalent
16 new points took8m37s, so allow roughly8–12minutes, without treating that
estimate as a completed measurement. Expect roughly another2.5GiB of artifacts;
the original936 packets occupy2,645,332,635 bytes and remain unchanged.

Whole-function reverse AST joins permit only the odd-index worker schedule,
the16/32 summary orders and the recorded old/new grid-address lookup. Every
full matrix, source binding, native quadrature, inverse, four-column solution
and end continuation body stays unchanged. Focused checks replay all16 accepted
addresses and full moment arrays exactly, reject wrong point addresses, and
cover every midpoint once. Reuse the accepted nonlinear-pole controls byte-for-byte.
No new integration is needed for these focused adapter checks.

Retain all375 analytic records,80 rows,35 sources,six profiles,160 native terms,
645 unknowns and four incident columns at every actual frequency. Keep frequency
dependence in complete end subspaces, incoming forcing, outgoing observation
and phases. Closed amplitudes remain boundary anchored. Physical row units and
the fixed positive seed scaling maps distinguish equation and boundary rows;
complex seed coordinates are not Hermitian flux-normalized physical currents.

Compare literal16/32 determinant winding increments and full inverse/response
moments0through3, all conditions and inverse residuals. Continue each complete
end cluster around the32-node circle and compare with radial node states and
loop closure. Preserve every point, denominator/gauge/Newton record and array.
Require clean workers/supervisor, empty stderr, checks/stdout identity, complete
native action/mass/measure/field/limit joins and all current/frozen/old/new hashes.
The16-point moments and increments must replay exactly from the unchanged packets.

This is the single selected refinement. Use its actual evidence to record a
scoped numerical outcome or pursue a candidate-local solve if a feature appears;
do not automatically launch another contour doubling. Numerical agreement does
not certify between-node holomorphy, all interior end exceptional loci, a global
empty spectrum or physical bound poles. A candidate requires full principal
parts and frequency-dependent overlap under nonlinearPoleV2, plus separate
sheet/width/normalizability/channel checks. No inverse moment is assumed a
projector, and zero residue cannot remove a higher-order pole.

On an instrument/runtime/schema problem preserve all original and completed
midpoint packets and partials. Resume only remaining work with exact helper and
source joins. Never repeat accepted integrations or solves for output plumbing.
One supervisor, durable repository scratch and a silent completion/error hook;
no model polling, routine updates, recurring checks or unrelated overlapping CAS.
