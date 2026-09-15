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

Implementation is committed at `185365fc`. The first constructor stopped
before saving a factorization packet: SymPy represents exponential powers as
`E: phase`, so looking for `exp(phase)` among power-dictionary keys rejected
a direct character. The repair reads the literal multiplicative factors and
uses a symbolic unit for a factor with no momentum character. Five focused
cases give 20 zero reconstruction/phase residuals, including repeated, mixed,
shifted and absent characters. The full native engine prefix has an exact
AST match to the accepted quadrature source. All original logs and snapshots
remain in the original run directory; no upstream reconstruction is needed.

The phase repair is committed at `91f56ec6`. Its retry processed the first
70 rows, then stopped at original index 70: the final ten native integrals
have middle momentum outside source position. Their native order is output,
input, source position, middle momentum. The other orders occur 40 times
(output, source) and 30 times (output, input, source). No complete or partial
factorization packet was saved by that run; its original sources and logs
are retained, and the upstream operand packets remain available.

The finite-domain adapter now identifies the unique source-position limit
by variable identity and preserves the relative order of every remaining
momentum limit. All 80 native layout checks pass. Six limit-order regressions
preserve the exactly evaluated bounded polynomial integral and replay saved
row objects; five phase cases also pass, for 32 zero residuals. The native
engine prefix remains an exact AST match. These are adapter checks, not
physical action convergence or a claim about infinite iterated integrals.

The next run preflights every layout and saves each completed row, phase
state and metadata context with hashes, before proceeding to the next row.
A strict source/provenance join permits resumption from those rows. It emits
the native source-limit index and remaining ordered limits with full metadata.
The limit-order repair is committed at `50c51500`. The all-operand run in
`retry-02` was launched with its owned supervisor and completion/error watcher;
no factorization output has yet been accepted or published. See the limit-repair checkpoint for the inventory,
regression records and original failure hashes.
