# Finish the one-row integration comparison

The principal-power continuation completed the full GK21 row in15.9346s with2667
quadrature evaluations, after one control point. It then hit the4096 cumulative
point cap inside GK15. All4096 full integrand input/value/completion packets,
source action, two controls and the complete GK21 return remain immutable.
Actual exits1 at26.529s native/26.703s guard; no interruption, cap/OOM/swap event.

Finish this same selected comparison with a bounded allowance of4096 additional
unique points (8192 cumulative), using the measured first-rule cost. Frequency,
finite domains, GK15 tolerances5e-9/5e-7, maxnorm and interval limit256 are unchanged.
This extends the point allowance for the unfinished comparison, not a frequency
or domain sweep. Keep the mandatory900s2GiBzero-swaponeCPU guard and measured cost
reserve. The original first-rule output and controls must never execute again.

The failed full quad_vec call has no final return or saved adaptive accumulator.
Its unfinished controller runs using the exact saved point function. Do not claim
that an absent adaptive state has been restored. Recombine those saved point
returns as needed by this unfinished comparison; do not reevaluate any completed
coefficient, frequency, Fourier source or full integrand. Reuse every prior point
only after exact variable/source-action/row/position input and byte-hash joins.
Preserve the new full result before its success guard. Record actual new and saved
point counts; do not invent the failed partial controller's unrecorded reuse count.

Whole original main, integrand and quadrature loop reverse-AST joins retain all
physical expressions, phase, source units, settings, principal branch and controls.
Only completed-prefix restoration, first-rule removal, exact existing input routing,
new metadata paths and the explicit numerical point cap differ. The final comparison
and full consumed input/source/logical/canonical hash guards are unchanged.

This remains one129x129 row at1-0.01i in an analog toy model. A bounded saved review
will assess the complete row comparison before a separate scoped pilot checkpoint.
Do not claim scattering accuracy, a global domain theorem or a pole result. Actual
four-case operators/four-incident solves and selected observable checks remain.
Preserve all failed and completed evidence; no automatic retry or old-directory
restart, no broad contour doubling or repeated ancestral science campaign.
