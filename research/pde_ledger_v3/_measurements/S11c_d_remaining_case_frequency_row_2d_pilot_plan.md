# One finite two-momentum row pilot

Use `exploratoryAcceptanceV1` for the analog toy model. The immediate question
is whether a compact numerical contraction can complete the first missing 2D
row economically enough to proceed to the full four-case scattering systems.
This stage computes only LAB_HELD__RHOBR_CONSTANT row46/source20, all129x129
entries, at the already chosen frequency1-0.01i. No new frequency sweep, root
selection, end map or global-domain certification is required.

The accepted two-row checkpoint is
`S11c_d_remaining_case_frequency_rows_1d_checkpoint.json`, SHA
2907d7da91fb048b345e0c4e4aa238a182b90e021672f81ecb453f7cf1fac2d5.
Its producer's176 memory.max contacts are retained; no OOM/swap occurred and
the hard2GiB cap held. The saved review was clean, including zero cap contacts.
Never rerun those numerical results to improve resource telemetry.

## Actual operands and new method

The saved row input contains one coefficient with literal Mul factors: a
constant group, an output-momentum/position group, an input-momentum group and
one actual finite nested profile Integral. Classify those existing factors by
their actual free symbols, without factor/expand/substitution, constructing a
new symbolic expression or omitting a term. Reject any other mixed factor.
The actual profile is
`Integral(10*(1/2-tanh(xi)**2/2)*exp(-10*I*xi*(k-q)),(xi,-14,14))`.
Require its literal `k-q` subtree and verify every occurrence of k/q in the
integrand is inside that subtree. Require the actual unit declaration and
retain the full original row/source/context/field-equation units/independent
grades/profiles/measure/Abel/settings/ordered limits and physical addresses.

Read the accepted saved prepared-source-basis.original arrays by complete
jet/probe/coefficient/unit/node/weight/size/caller identity. The source action
has1024 rows and129 columns; it already includes the source measure. Route
matching actual row29 Fourier inputs/values/receipts, with original
computational owner separate from the new consumer. An absent serialized
internal return is not permission to recreate a completed old call. No old
source_jets, derivatives, symbols, Poly, compiler, lambdify, prepare_basis,
Gauss-rule constructor, current, mode, end map or baseline calculation runs.

Compile only the previously reviewed numerical expression-tree evaluator,
with the accepted native principal-complex power convention. This is new
numerical evaluation of exact saved expressions, not symbolic derivative or
compiler reconstruction. Preserve native pickle restoration; do not patch
Basic class methods. Unsupported nodes stop with inputs already saved.

Use a new explicitly labelled uniform composite trapezoid on[-4,4]^2 with512
and1024 panels per leg. This differs from the original16/4/4 split Gaussian
rule, which remains recorded. Dyadic spacing gives exact repeated difference
addresses. Evaluate the finite profile integral on a new4096-panel uniform
trapezoid over[-14,14], in batches of at most64 differences. First compare
five actual differences0,+/-.731,+/-8 against a new2048-panel profile rule,
absolute target1e-9. These are new rules, not claimed restoration or completion
of native profileOrder512. Physical cutoffs/regulator remain unchanged.

For each genuinely missing q compute the literal minus-phase Fourier vector
against the saved weighted source action. Reuse exact completed vectors for
matching full calls and nested grids. Evaluate each literal coefficient factor
only on new grid coordinates. Constants and exact coarse values are reused.
With A(k,z), B(q), C the actual numerical factor products and P(k-q) the saved
profile values, calculate the complete row as
`A.T @ (wk[:,None]*P*wq[None,:]*B[None,:]*C) @ Phi`.
This is a regrouping of all finite quadrature terms; T is an ordinary transpose,
with no conjugation. Preserve input/intermediate/value/completion records
before subsequent guards, including the profile integrand, kernel, inner
contraction and final full row. Use compact array packets rather than storing
a129x129 outer product for every momentum pair.

## Focused checks and stopping

Before contractions require new rule mass/orientation, finite complete arrays
and actual native source/unit joins. Compare three actual full coefficient
evaluations against the grouped factors at off-diagonal coarse-grid pairs,
scaled target2e-12. Compare nine selected row entries with literal double sums
from the same saved arrays, scaled target2e-12. Require actual changed-measure
and transposed-kernel controls to respond. These new algorithm checks do not
repeat an old evaluator, source, quadrature or completed proof.

Report full-row maxnorm and coarse/fine absolute and relative spread. The pilot
resolution indication is absolute1e-4 or1% of row norm. This is not an
observable error bound or evidence resolving a tiny physical effect. Preserve
an unresolved indication if this comparison misses the target; no automatic
doubling, alternate tolerance or campaign follows. The eventual scattering
observables still need selected comparisons near1%, amplitude1e-4/current1e-6.

The previous 1D integration costs were24.46--46.47s each. The proposed compact
2D route has about10million finite profile samples and at most1025 source
Fourier frequencies, with matrix contractions in bounded arrays. Its cost is
initially unmeasured, not inferred from source-family counts. Measure the first
complete512-panel row; run the1024-panel comparison only if3x that measured
cost plus60s fits the remaining900s. Record actual costs and all outcomes.

## Containment, publication and continuation

Use one guarded worker:900s whole-job limit,2GiB,zero swap,one CPU,nice15,
32tasks,nativeThreads1,no OOM. No overlap, automatic retry, old-directory
restart, unguarded fallback or model polling. Keep a silent local completion/
error hook. New packet writes flush their own bytes before advising away clean
pages; streaming reads also advise clean pages away. This changes no file
bytes or resource cap and makes no claim about the old cap-event cause.

Require final actual guard/supervisor/child0, empty strict stderr, stdout/checks
byte identity and consumed current/frozen/source/logical/canonical/hash/size
postchecks. Preserve all failures and every completed input/value/receipt.
Bounded saved review must precede scoped acceptance, without recomputation.
Then use actual cost/resolution to choose remaining row work toward645x645,
four-incident, four-case frequency systems and exports. No whole-program,
physical-pole or outgoing-domain completion is claimed.
