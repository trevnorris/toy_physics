# S11c-d source Fourier quadrature checkpoint

The bounded factorization is published and annex-verified at bed088be: all 80
native operators and 35 source integrals, with 324 zero normalized residuals
and 264 zero proof residuals. Original denominator/branch restrictions and raw
representation residuals remain explicit.

The source quadrature adapter binds every source factor to both approved
Gaussian tests, derives its actual frequency range, and compares three Gauss
orders with literal original-source and independent adaptive integration at
two finite source intervals. Phase batches have an explicit memory budget.
Each binding and completed integral is saved for recovery. Focused finite
integration errors are below 3.56e-15 and are unchanged by batch size.

The native binding comparison now follows all 320 field/row/term occurrences.
Ten original/test pairs differed live but all saved occurrences compare exactly
equal after pickle reconstruction; live representation strings/hashes remain
alongside raw residuals and exact certificates. The first loader failure and
the later exact-tree join failure stopped before numerical quadrature. Their
logs and operands are preserved in the execution and repair checkpoints.

An expanded-form CSE stress check entered an eight-hour polynomial GCD. The
bounded uncompressed alternative passed all ten cases in 75.22 seconds, with
peak worker RSS 71,708 KiB: all exact residuals and 62 proof scalars vanish,
all native limits agree, and every coefficient mutation remains nonzero.
All saved operand/source/hash joins passed validation. Only the production
comparison's shared keyword changed to False; the complete checker AST joins
after that substitution and the physics engine is unchanged. The repair is
committed at 420f36a0.

Full source quadrature retry-02 is running with the single-job supervisor and
silent local completion/error watcher. No full source-quadrature result is yet
accepted. Next is validation/publication, followed by concentration-aware
full-action momentum integration. Finite source tests establish no uniform
interpolation bound, infinite tail/interchange, full-action convergence, Abel
weak limit, scattering or pole result. See the plan, execution checkpoint and
binding acceptance record for evidence and retained boundaries.
