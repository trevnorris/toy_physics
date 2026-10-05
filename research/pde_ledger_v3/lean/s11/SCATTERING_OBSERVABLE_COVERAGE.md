# Scattering observable work contract

Authorized 2026-09-26. Governed by ../FORMALIZATION_POLICY.md. Status: bounded F1–F4 complete after run8 verification and independent
Claude/Grok CLEAR reviews. See SCATTERING_FLUX_FIDELITY_REVIEW.md. No production calculation is authorized
by this document. The current native response/current construction is evidence to
identify conventions, not a premise that its unfinished physical boundary problem
has been solved.

## First increment: F1–F4, finite current and flux algebra

| Item | Bounded claim and coverage |
|---|---|
| F1 | On every finite complex amplitude space (including dimension zero), use the entire supplied matrix J: pair(J,x,y)=conj(x)^T J y. Prove its exact coordinate expansion, conjugate symmetry when J is Hermitian, and the full quadratic interference identity. The real flux is Re(pair(J,x,x)); no arbitrary imaginary physical-current defect is certified away. |
| F2 | Prove pullback by C is C^H J C, preserving all entries, and scattering-coordinate covariance under supplied inverse basis maps. Rectangular pullback alone is not a basis equivalence or a complete channel census. |
| F3 | Retain left/right outward signs -1/+1 and the opposite incident sign. A normalized fraction is undefined at zero incident flux; negative, zero and positive flux exhaust the real domain. Nonnegativity and a bound by one require explicit numerator/denominator/balance premises. Sector and incident/outgoing cross terms disappear only under an explicit vanishing premise. |
| F4 | Exact nonzero examples and paired true/false controls test conjugation, cross terms, basis normalization, orientation, zero denominator and conditional balance. No invalid import, warning, timeout or resource failure counts as rejection. |

Fidelity: directives/S11c_d_SHARED_PHYSICS.md sections 3a/3c and selected
functions in _measurements/S11c_d_continuum_currents.py and
S11c_d_continuum_response.py. Only small selected algebraic functions may execute
in an isolated native check; no module import, saved scientific packet load,
constructor or production rerun. Source hashes and explicit normalization are
recorded. The native translation remains outside the Lean kernel. Physical flux
units are supplied by J; ratios require numerator/denominator in the same units.
No factors of frequency, epsilon or Fourier measures are inferred from this algebra.

Completion: canonical proof/audit, compact native identification, mathematical
controls and positive examples, current source/dependency/object correspondence,
two independent non-author fidelity reviews with findings resolved. Stop there.
External packet transfer requires authorization for the new packet; prior
increment approvals do not cover it. No commit requested.

## Local evidence and case coverage

| Obligation | Proved contract | Controls and evidence |
|---|---|---|
| F1 | `pair_coordinates`, `pair_hermitian`, `hermitian_flux_real`, `flux_add`, `flux_add_iff_cross_zero`; every finite complex space, all matrix entries | Interference and conjugation pairs; non-Hermitian imaginary witness and false reality claim; empty-space and nonzero indefinite-null positives |
| F2 | `pair_pullback`, `flux_pullback`, `pullback_comp`, `scattering_flux_covariant` with its supplied right inverse | Basis-metric pair; native nonunitary complex congruence, wrong transpose and stale-metric checks |
| F3 | `end_coverage` and `flux_sign_coverage`; `fraction_undefined_iff`/`fraction_defined`; conditional bounds and balance | Orientation, signed incident, zero denominator, normalization, negative denominator, false unconditional upper bound and false unconditional conservation pairs; balance/zero-defect positives |
| F4 | Exact examples in `S11ScatteringFlux/Controls.lean`; fresh true/false statements decided using canonical witnesses in the instrument | Eleven mathematical rejections, fifteen positives, each rejected statement reaching exactly one literal `False` in `contract_control`; no source-replacement mutations claimed |

Left/right are the two distinct constructors of `End`. The strict negative,
zero and strict positive alternatives are mutually exclusive by the real order
laws; each is inhabited by the supplied witnesses. Zero dimension is permitted.
Zero flux includes nonzero amplitudes of indefinite forms. For nonzero incident
flux, negative and positive denominators remain separate physical-interpretation
cases; the algebra alone imposes no positivity. Cross-term cancellation is the
zero versus nonzero real sum of both contractions, not a diagonal-mode assumption.

Recorded run8 freshly built three modules and the audit root, checked 34 selected
standard-axiom lists and passed all 26 controls. The compact native report passed
seven identities/examples and five wrong-formula controls. See
[SCATTERING_FLUX_VERIFICATION.txt](SCATTERING_FLUX_VERIFICATION.txt) and
[SCATTERING_FLUX_FIDELITY.md](SCATTERING_FLUX_FIDELITY.md). Local execution and
source/object correspondence was revalidated at closure; both independent
source-fidelity reviews are complete. Neither reviewer reran the Lean build.

Excluded: physical mode existence, completeness of any supplied numerical channel
basis, reality/positivity of actual derived currents, an actual scattering solve,
unitarity or conservation without its balance premise, radiation/centre boundary
repairs, bound-state capture, infinite-dimensional convergence and CAS bridging.

## Authorized subsequent increments

1. Retained-order bookkeeping: exact independent eta/sigma rectangle and its
   named one-parameter path; quadratic flux coefficients through order two,
   changing current forms and incident denominator, total versus coherently
   subtracted induced flux. Formal coefficients are not parent-theory Taylor
   coefficients without a remainder premise. Write the bounded follow-on
   contract before implementation; use F1–F4 rather than re-proving its algebra.
2. Finite-solve sensitivity: residual-to-amplitude-to-flux estimates using the
   existing analytic-error stability lemmas, once the physical finite solve and
   norm/observation maps are identified. A measured condition number is not a
   certified inverse bound. This is not a gate to exploratory numerical work.
3. Portable D5 bulk registration: add the already independently reviewed D5B
   target, validate import/control/native-input census and tooling regressions.
   Preserve historical INSTALL validation and proof records. Do not repeat the
   expensive accepted D5 proof solely to test a catalog entry; distinguish
   registration validation from a new full fresh-object execution.

Execution uses the host systemd resource guard, one job at a time, at most
2 GiB whole-job memory, no swap, one CPU and low priority. Lean probes use one
worker and a smaller allocator bound. Failure of containment stops execution;
no unguarded fallback or automatic retry. Preserve all historical and S11c files.
