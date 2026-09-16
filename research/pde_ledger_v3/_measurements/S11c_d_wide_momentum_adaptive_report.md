# S11c-d independent wider-box outer quadrature

Production completed cleanly in 4m20s with four single-thread workers.
Independent GK21 outer quadrature agrees with the refined Gauss reference
within 4.83e-15 for all40 single-momentum rows and 7.71e-16 for all30 paired
rows on both fields and cutoffs3/4. Complete-action differences reach
4.76e-15 / 3.47e-18. The cutoff3-to4 action change remains 1.134674e-7.
All eight solves finish with reported unit-frame error estimates below9.58e-12
at requested absolute tolerance1e-10 and zero relative tolerance.

Position bounds48/14, regulator0.2, source/profile256/512 and paired inner64
were held. All ten three-momentum rows remained explicitly at144/24/24;
Gauss outer432 was only the reference for the independently integrated layouts.
This rules out observed single/pair outer quadrature error as the source of the
much larger sampled box difference. Wider-box three-momentum refinement is next.

Saved-operand acceptance verifies all84 current/frozen sources, 3066 conditional
points, 3090 worker artifacts, 12 records, 12 exact read-only source/profile
caches and2160 zero held-term scalars. All native row/source/profile/limit and
conditional-unit joins, actual conditional-weight mutations, finite masses,
completed adaptive intervals, literal term/action differences and packet hashes
pass. Full replay covers14882 tags,7439 fresh keys and101909 metadata paths.
Workers/supervisor exited zero with empty stderr and checks/stdout identity.
Peak worker RSS was131488KiB. No integration was repeated for acceptance.
Preflight5d519d3a remains instrument evidence only.

The canonical6,750,781-byte transcript is
`scripts/out/S11c_d_wide_momentum_adaptive.out`, SHA256
`afdac193d8b7d3fce4bb734336503d1d02d182e38f40829d914261a97d4be99d`.
Full wide-box/independent-grade/uniform/global exceptional convergence,
physical infinite tails, Abel limits, scattering and poles remain open.
The physics engine, approved values and pinned exports are unchanged.
