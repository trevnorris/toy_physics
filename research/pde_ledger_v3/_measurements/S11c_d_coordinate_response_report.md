# Affine material scattering route

Coordinate sources are accepted and annex-verified at9c4382c6. The numerical
adapter now binds actual accepted coordinate images, the complete transported
field-derivative basis, all80 integral rows and35 source amplitudes. It maps
native momentum nodes and measures, source/profile limits and their Jacobians.
All independent local/cell coefficients remain separate. No physical input or
native engine changes, and no existing S-matrix is conjugated.

The bounded focused run compares actual material and Eulerian row integrands
on production-order source/profile rules and small momentum prefixes. It also
constructs the material end trace/insertion/phase maps and transforms the52
saved polarized-current derivative tables before full-subspace contractions.
Data return to the common Eulerian coordinates inside construction. Actual
missing-Jacobian and wrong-phase controls must respond. This first run has
not yet established numerical agreement; its results must pass final guards.

The complete material momentum quadrature and coefficient response are the
next implementation step after the focused check. Reuse saved binding and end
operands with exact source/helper joins. The full material mass guard must use
the actual transformed box and Jacobian; it must not reuse the unmodified
physical-box volume check. Retain approximate boundaries, finite regulator,
transported trial-space scope and practical observable resolutions. The
one-sided shape scattering sensitivity, targeted frequency poles, remaining
cases and final exports remain. No broad convergence campaign.
