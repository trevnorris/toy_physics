# S11c-d local and nonlocal assembly

Validated on the supplied LAB_HELD / RHO4_CONSTANT source. The constructor
computes four local derivative matrices (orders 0–3) and retains 80 distinct
nonlocal integral operators across all 25 row/column actions: 71 local and
160 nonlocal coefficient entries. Every original output/input/middle
momentum argument, profile transform and integration limit remains on its
operand.

All 281 scalar residuals are zero: full-action reconstruction, affine/nonlinear
remainders and coefficient differences against an independent derivative
extraction in the formal carrier algebra. The 170-tag transcript passed full
saved-object replay, with 1,226 resolved dimension/grade paths and 83 fresh
injective write-keys. The constructor is ReducedActionAssembly.construct in
the existing engine; S11c_d_reduced_action_assembly_check.py performs the source,
cache, emission and residual checks. These establish exact assembly, not
numerical quadrature, convergence or scattering.

The initial run saved its packet and transcript, then failed at JSON summary
serialization of a SymPy derivative-order integer. Commit 70d605c2 fixes that
instrument issue. Recovery reused both files byte-for-byte, verified the
original source snapshots and unchanged loader/emitter, and completed final
validation in 50.61 seconds with empty stderr (154,500 KiB peak RSS). No physical
formula or assembly calculation was rerun. Both original and recovery logs and
provenance remain in repository _scratch/s11c/.

The 829,611-byte accepted transcript is published as
scripts/out/S11c_d_reduced_action_assembly.out through DataLad/git-annex.
Its SHA-256 is fe1d8d656621298e056ea4078312d353fc0e0b4a2742f03636413b5c38a47658.
The complete evidence inventory is S11c_d_reduced_action_assembly_checkpoint.json.

Numerical action evaluation needs values for 30 inherited free gradient-energy
coefficients missing from the endpoint-only input. S11c-b section 3a explicitly
leaves these constants free; this is an input choice, not an upstream bug.
Their symbolic dependence is preserved. Await the requested numerical-instance
choice before binding them. Then verify the extended input's end restrictions
and compute actual local/nonlocal actions against direct reduced-row test
actions, with quadrature, tails and the Abel weak limit accounted for before
boundary matching. The S-matrix, continuum expansion, profile-frequency poles
and final own-row export remain open program work.
