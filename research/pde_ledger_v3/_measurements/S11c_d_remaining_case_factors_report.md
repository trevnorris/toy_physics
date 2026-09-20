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

The first recovery stopped before certification because its generated global
helper name collided with the worker local `certificate` variable. A three-site
alias correction has a reverse whole-file AST join and two verified global
call sites. The certificate algorithm and all three completed rows are unchanged;
the fourth retained raw row remains byte-identical. Full execution/output is
still pending, with the failed recovery logs and frozen sources preserved.

All four workers now pass. Row 28 completed its saved-row certificates in 4.63
seconds with two certificates, ten zero proof scalars, four zero character
scalars and a responding mutation. The supervisor then rejected one metadata
path: provenance unionIndex=0 lacked the dimensionless zero-unit annotation.
All factor and proof packets are complete and remain unchanged.

The output adapter supplies only that association path's dimensionless unit.
The actual saved payload is unchanged; three changed-unit controls and a wrong
path control reject. Its coordinator reverses exactly to the prior implementation
after four output/provenance/resume wiring changes. The native factor, certificate
and emitter/checker sources are unchanged. Full saved-row output acceptance is
pending; no factor, certificate or integration will be recomputed for this repair.
