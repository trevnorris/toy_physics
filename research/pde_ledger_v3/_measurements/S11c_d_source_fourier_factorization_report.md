# S11c-d bounded source Fourier checkpoint

The quadrature-domain result is published and annex-verified at `e514b704`.
The next numerical assembly stage separates source-position integrals from
the remaining momentum quadratures, using the accepted original nonlocal
operators and preserving their full symbolic parameter dependence.

BoundedSourceFourierAssembly in the existing engine introduces explicit finite
cutoffs and collects source-dependent factors, source-independent coefficients
and their actual exponential characters. It derives the source frequencies
from those characters and constructs bounded source Fourier operands. The
checker retains all original integral/field/branch joins and tests literal
integrand, amplitude and phase residuals with full dimension/grade replay.

This is a finite-domain construction for the stated smooth-test, positive-
regulator and nonsingular-denominator conditions. It does not exchange the
original infinite distributional integrals or establish their physical limit.
No source integral is numerically evaluated by this factorization alone.
Validation is pending before applying the split momentum panels and advancing
the full-action convergence checks.

Implementation is committed at `185365fc`. The native one-case constructor/
validator is running with its owned supervisor and silent completion/error
watcher; the factorization result is not yet accepted.
