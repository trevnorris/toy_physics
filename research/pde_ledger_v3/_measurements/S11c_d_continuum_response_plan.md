# Complete finite continuum coefficient response

Interior matrices are accepted at7d1082df and actual boundary/current
coefficients at894420f6. Combine their common approved input and finite
interval/source64, momentum4, profile14, regulator0.1,129-coefficient basis.
No new quadrature or end-mode construction is needed.

Replace only the actual endpoint rows of every interior coefficient matrix
with the computed trace-map coefficients. Insert the complete incoming boundary
data for both ends. Solve the baseline and independent eta,sigma,mixed equations
recursively from the same baseline operator. The mixed forcing includes both
cross terms and the actual mixed operator, with no frozen boundary substitution.
Save all four systems, right-hand sides, forcings, solutions and residuals.
Use equilibrated LU and an independent SVD route, full rank/condition checks,
and a one-sided mixed-forcing omission control.

Extract every outgoing modal amplitude, including the three evanescent matching
directions at each boundary. Restore a common origin for open amplitudes using
the computed phase coefficient maps. Keep evanescent amplitudes at their real
boundaries. Provide both the reference-anchored gauge and the field-coordinate
map derived from the accepted reference normalization, with source classifier
labels. The field gauge is analytic continuation of the reference field basis.

Use the complete computed open-channel current coefficients to construct their
positive square roots by coefficient Sylvester equations. Derive the current-
normalized S-matrix from those maps, including their eta/sigma dependence.
Keep baseline/interference/quadratic homotopy coefficients of open outgoing
current, its incident-current denominator and their quotient. These describe
the open-channel response; do not call an evanescent matching amplitude an
outgoing thickness flux, or infer bulk escape from an open-current deficit.
Thickness/bulk conversion, physical controls and poles remain separate work.

Check every coefficient equation and endpoint trace, independent solve,
normalization square/inverse/Hermitian residual, actual source/field/channel
join and phase map. Compare the coefficient solution with direct solves of
the same assembled coefficient polynomial at selected small formal points;
record the retained-expansion remainder, not a fit used to construct coefficients
or a physical parameter change. Do not identify that remainder with an error
bound for the original parent theory.

Emit compact full-array fingerprints and literal per-field residuals with units,
independent grades and lambda order. Preserve full coefficient arrays and all
source hashes. Replay every payload and metadata path. One native thread,2GiB
and900seconds bound the attempt. Publish and commit accepted evidence, then
finish physical conversion/current bookkeeping, controls, targeted poles and
remaining cases/exports under the practical toy-model standard.
