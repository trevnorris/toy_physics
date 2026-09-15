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

The first launch stopped before binding or numerical integration: the source
checkpoint contains two publications, while the new loader expected one. The
loader now validates both existing schemas and every listed publication. All
four consumed checkpoint inventories pass; only that loader changed. Original
logs and source snapshots are preserved. The retry stopped before quadrature on the exact native-test binding guard.

The binding diagnostic covers 160 original/test pairs and all 320 explicit
field/row/term occurrences. Ten pairs had unequal live expression trees; all
limits agree. Pickle reconstruction makes all 320 saved occurrences exactly
equal. Certificates on those restored pairs have zero normalized residuals
and 52 zero proof scalars; a coefficient mutation remains nonzero. The repaired
production helper also rejects a changed integration limit.

The checker now compares each actual occurrence, preserves raw and certified
residuals plus live representation strings/hashes, and saves per-test joins
and each bound source before later guards. It uses the existing exact
certificate helper; the engine and physical inputs are unchanged. The expanded-form regression ran for more than eight hours without a completed
result. SIGINT preserved a traceback in multivariate polynomial GCD inside
`sp.cancel`; SIGTERM then finished interpreter teardown. No full source quadrature has
yet run or been accepted. See the binding repair checkpoint for saved operands,
source hashes and the live regression record.

A bounded replacement tests the existing uncompressed (`shared=False`)
certificate path. Each of the ten saved cases runs in its own subprocess,
with a three-minute wall limit, an 8 GiB address-space ceiling, an early
traceback, and operands/certificates saved before subsequent work. Every case
must pass along with its coefficient mutation; a timeout remains unresolved.
The production certificate criterion and physics engine are unchanged pending
this test. The accepted source/factorization/test packets remain intact.

The bounded uncompressed regression is validated: all ten cases completed in
75.22 seconds, with peak worker RSS 71,708 KiB. Every expanded representation
differs from its comparison operand, all ten exact residuals and all 62 proof
scalars vanish, all original limits agree, and every coefficient mutation is
nonzero. Saved certificate/mutation operands, case/source hashes, and all
stdout/checkpoint joins passed inspection.

Production now selects `shared=False` only in `integral_comparison`. A whole
checker AST join permits exactly that one boolean change; the physics engine
is unchanged. Raw residuals, live representation strings/hashes and all 320
occurrence joins remain. Full source quadrature is prepared for retry-02.
