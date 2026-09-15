# S11c-d two-momentum refinement

The one-momentum refinement is published at 9fad1d11. Its final outer/source
changes are below 6.64e-15/3.81e-15, and adaptive outer actions agree within
1.90e-15. This stage retains those values and the original finest three-momentum
terms while refining all 30 native two-momentum rows on both approved fields.

The prepared run varies outer order, inner panel order, source order and profile
order separately, then compares independent adaptive outer quadrature. It keeps
the original finite bounds and positive regulator. Complete actions explicitly
identify the changed and held rules; held-layout convergence is not presumed.

The new helper caches the same source and profile Gauss nodes/weights. A whole-
engine AST join preserves all existing construction code, and a separate method
join permits only the profile-rule lookup change. Every cache is compared exactly
with the original rule and made read-only. One native numerical thread is used.
Every completed grid, partial sum and adaptive outer-point value is saved.

Focused checks passed in 31.9 seconds with exit zero and empty stderr: all
30 rows and complete actions reproduce the accepted baseline exactly. Summing
128 independently evaluated inner-at-outer values differs by only 1.17e-17.
All source/profile node and weight caches match exactly and remain read-only;
actual weight mutations respond on both tests.
No paired refinement result is accepted yet. Three-momentum refinement,
independent-grade coverage, physical tails/interchange and Abel weak limits
remain separate requirements before matching, scattering and poles. All approved
inputs, durable operands and the retained solver/export contract are preserved.
