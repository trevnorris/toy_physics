# S11c-d source Fourier quadrature checkpoint

The bounded factorization is published and annex-verified at bed088be: all 80
native integral operands and 35 source integrals are retained, with 324 zero
normalized residuals and 264 zero proof residuals. Raw representation
residuals and original denominator/branch restrictions remain explicit.

`BoundedSourceFourierQuadrature` evaluates actual bound source amplitudes in
finite intervals with bounded phase batches. The checker consumes the
accepted source/certificate/test packets, derives frequency ranges from each
actual affine phase, and compares all source transforms on both approved tests
with the original bound source and independent adaptive integration. Three
Gauss orders and two finite source intervals are retained separately. Each
completed integral is saved, hashed and resumable with source/provenance joins.

Focused finite-integral errors are 3.56e-15 (Gauss) and 1.78e-15 (adaptive).
Changing the batch size changes no result; the native engine AST is unchanged
apart from the new class. The full accepted-source run is pending. No new
physical input, upstream repair, full-action convergence, infinite-domain
interchange, Abel weak limit, scattering or pole result is established here.

Next: validate/publish the source quadrature and use it in concentration-aware
full-action momentum integration. See the plan and execution checkpoint.
