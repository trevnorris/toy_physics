# S11c-d quadrature domain checkpoint

The source-derived resolution check is validated for the approved LAB_HELD /
RHO4_CONSTANT instance. It retains 23 coordinate-dependent denominator bases
and all three Abel momentum-transfer pairs. Their computed real-part sign
records retain the positive-regulator/real-momentum assumptions. An unset
SymPy sign flag is not the opposite property or a global spectral certificate.

The existing Abel kernel gives a physical momentum half-width a/10. Its
computed primitive and half-height residuals are zero. Ordinary Gauss rules
remain poorly resolved at narrow peaks; across the sampled regulators and
centers, the maximum bounded-mass error is 0.03768 even at order 1024.
Rules split around the computed width reduce that error to 5.33e-15 with
order 16 on each panel. All centers, regulators, panels and operands are saved.

All six native profile transforms retain their original momentum-leg maps.
Split Gauss integration at order 128 per half-line agrees with adaptive
finite integration to 1.25e-12 across the sampled transfers and cutoffs.
At profile cutoff 10, estimated two-sided tails are 2.07e-9 for the zero-jet
subtraction and 4.13e-8 for the derivative transforms; at cutoff 14 they are
6.92e-13 and 1.39e-11. These are elementary profile estimates, not bounds on
the complete operator action or an established Abel weak limit.

The first run preserved all operands but cast its complex Abel sums to real
numbers, producing two warnings. Repair 84968950 retains full complex values.
Recovery reused every denominator/profile construction and recomputed only
small sums from the saved density and quadrature operands. Every old real
projection agrees exactly; restored imaginary parts reach 1.58e-16. The
original denominator packet is byte-identical, profile objects are reused,
and all seven consumed-helper joins pass. Original artifacts remain intact.

Recovery completed with exit zero and empty stderr in 45.81 seconds. Full
replay covers 310 tags, 6,222 metadata paths and 153 fresh write-keys; the
saved packet hash is unchanged across emission. Source, artifact and stdout
joins were verified again before publication. The 307,234-byte transcript is
published as scripts/out/S11c_d_quadrature_domain.out through DataLad/git-annex,
with SHA-256 e31fc2bb175182198135b3f512935dd37a34ab4193878a79382f56468fe6f014.

Next: apply the source-derived panel scales to full-action refinement with
bounded memory. Separate the source-coordinate Fourier integrals on finite
domains before constructing larger momentum quadratures, retaining literal
factorization/phase residuals and original leg/measure provenance. Full action,
domain, tail and Abel weak-limit convergence still precede boundary matching.
