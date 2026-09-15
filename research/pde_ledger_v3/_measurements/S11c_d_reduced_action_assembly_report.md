# S11c-d local and nonlocal assembly

The reduced-action source is committed and annex-verified at `83e98566`.
`ReducedActionAssembly` now collects the actual native probe columns into
local derivative coefficient matrices and intact field-dependent integrals.
Output/input/middle momentum labels, inner profile transforms and every
integration limit remain on their source operands. No integral is evaluated by
polynomial collection or by its fingerprint.

For each of the 25 actions, the constructor computes a literal reconstruction
residual, any affine/nonlinear remainder, and coefficient residuals against an
independent differentiation of the source in the formal carrier algebra.
The checker preserves the full packet before emission, verifies all source and
cache joins, and replays every fingerprint, dimension and grade. Results are
pending the run. These checks concern exact assembly, not physical quadrature.

The source contains 30 inherited free gradient-energy coefficients absent from
the current numerical input. S11c-b section 3a explicitly carries these constants;
this is not an upstream implementation repair. Their restored units and source
components are inventoried in S11c_d_variable_profile_parameter_inventory.json.
A numerical-instance choice has been requested. Keep all coefficients symbolic
until that choice; symbolic assembly does not depend on their numerical values.

Next: bind the complete agreed parameter map and independent profiles, compute
local and nonlocal actions against direct reduced-row test actions, and examine
quadrature, tails and the Abel weak limit before boundary matching. The full
S-matrix, continuum expansion, profile-frequency poles and final export remain
open program work.
