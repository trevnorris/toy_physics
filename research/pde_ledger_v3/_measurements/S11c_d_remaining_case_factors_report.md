# Remaining case Fourier factors

Accepted all-case source comparisons identify 32 distinct new integral operands:
9 first occur in LAB_HELD/RHOBR_CONSTANT, 19 in MATERIAL_ADVECTED/RHO4_CONSTANT,
and 4 in MATERIAL_ADVECTED/RHOBR_CONSTANT. The last case also uses 20 new
operands already present in the other cases. The native factorization and its
transitive helpers join the accepted frozen source exactly. No new operand has
three momentum variables.

The preflight saved all four native factorizations. Three rows completed cleanly
in 2.28, 3.21 and 0.86 seconds. The fourth worker stalled in polynomial GCD
while checking an unreachable Piecewise branch, then exited -11; its complete
raw row and live pair remain saved. No unfinished preflight is accepted.

The local certificate repair proves each branch's first-match condition and
omits only an impossible effective condition. It uses exact rational numerator
expansion for the reachable identities, retaining original denominators and all
raw branches. Four saved-pair checks pass in 6.30 seconds with 26 zero proof
scalars and four responding coefficient mutations. Changing the actual earlier
guard exposes the nonzero shadowed branch and is rejected. Native encoder and
worker changes have reverse whole-function AST joins; the engine and native
factor constructor remain unchanged.

Recovery reuses the three completed row/proof packets and the fourth raw native
row, computes only its missing certificates and then finishes output/metadata
validation. No source reduction, factor construction, numerical integration,
mode or solve repeats. The remaining 28 unique factors will follow only after
this preflight passes. Remaining case responses and final exports are still work.
