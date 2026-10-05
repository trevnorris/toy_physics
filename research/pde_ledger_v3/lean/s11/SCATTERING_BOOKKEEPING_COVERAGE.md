# Retained-order observable bookkeeping: P1–P4

This is the second increment authorized on 2026-09-26 in
SCATTERING_OBSERVABLE_COVERAGE.md. Status: bounded P1–P4 complete after recorded
run7 verification, author correspondence validation and independent Claude/Grok
CLEAR reviews with no blocking findings. Work is paused at the user's request.
F1–F4 is already complete with both reviews CLEAR.
Follow ../FORMALIZATION_POLICY.md. No new physical solution is assumed.

| Item | Claim |
|---|---|
| P1 | Keep eta and sigma independent in the four-term amplitude rectangle. Under the supplied real path eta=lambda, sigma=r lambda, prove a0=a00, a1=a10+r a01, a2=r a11. Delta includes both zero-jet contrast and first-jet pieces. The path does not identify the two formal grades. |
| P2 | For actual finite amplitude polynomial a0+lambda a1+lambda^2 a2 and supplied current polynomial B0+lambda B1+lambda^2 B2, prove the coefficients of the full quadratic contraction. Baseline, both coherent interference terms, current variations and a2 terms all remain. An exact higher-order remainder records terms above degree two rather than discarding them by notation. |
| P3 | Given nonzero incident coefficient j0, prove the quotient coefficient equations c0=n0/j0, c1=(n1-j1 c0)/j0 and c2=(n2-j1 c1-j2 c0)/j0. State separately the zero-j0 case. These are coefficients of the supplied retained rational model; no parent-theory remainder bound or nonzero denominator on an unspecified interval is asserted. |
| P4 | Compare the total contraction, total minus baseline and the coherently subtracted induced quadratic form. Prove the baseline-free specialization and exhibit nonzero-baseline counterexamples. Mathematical controls target omitted current variation/interference/a2/incident-denominator terms, wrong path ratio and conflated induced/total observables. |

Fidelity anchors: shared physics sections 3c/3d and the existing selected native
`amplitude_bookkeeping`, `lambda_series`, `quadratic` and `quotient` functions.
Their actual finite polynomial/convolution operations are the compact object
link. No physical packet, solver, production module or multidimensional sweep
needs to execute. Reuse F1–F4's conjugation and full-current conventions. Explicit
real paths, finite complex amplitudes and all coefficient signs are retained.
Any interpretation as physical real flux requires the corresponding reality
premises; any interpretation as a fraction uses the incident-flux domain.

The epsilon amplitude factor gives epsilon-squared flux; cancellation in a
fraction requires nonzero epsilon and nonzero incident flux. The retained
rectangle does not supply omitted pure eta^2/sigma^2 parent-theory amplitudes.
With nonzero baseline they can change the parent-theory second-order observable;
this increment must not promote a retained coefficient to that stronger claim.

Completion: proofs of the declared finite algebra, explicit domain/case coverage,
compact native identification, positive and mathematical negative controls,
standard-axiom audit, source/object validation and two independent non-author
fidelity reviews. Excluded: full parent-theory Taylor/asymptotic claims, strong
contrast extrapolation, actual conversion values, physical boundary repairs,
bound capture, scattering convergence and CAS transcript replication.

## Concrete implementation and fidelity limits

`Path`, `Quadratic`, `Quotient`, `Observables` and `Controls` implement the bounded
contract. P2 uses the complete supplied degree-two amplitude/current polynomials;
its exact complex contraction has degrees zero through six, and the real-flux
statement follows by taking real parts. Empty finite amplitude spaces are allowed.
P3 is field algebra, applicable both to real flux coefficients and to the native
complex diagonal channel coefficients. A complex quotient is not automatically
a physical real fraction. Each supplied incident column is treated separately;
all off-diagonal current entries inside its quadratic numerator remain.

The native `lambda_series` maps every supplied bidegree (p,q) to degree p+q
with weight r^q. `current.multiply` is unrestricted convolution; the response
engine's componentwise cutoff is not used for this native product check.
Selected original `amplitude_bookkeeping`, `subtract`, `quadratic`, `multiply`,
`quotient`, `lambda_series` and `adjoint` functions may execute only on small
synthetic fixtures. No production import, saved scientific operand or driver runs.
Record Python/NumPy versions as well as source/AST hashes.

The original native quotient has no zero-leading-denominator guard and preserves
complex values. Its correspondence checks use the declared nonzero domain;
an isolated zero-denominator fixture produced nonfinite output, demonstrating
that an external domain check is required. This does not establish that any
saved physical denominator vanishes.
When j0=0,n0≠0 there is no ordinary power-series constant solution. When both
are zero, existence/uniqueness must be resolved separately; the all-zero example
has nonunique coefficients. This contract supplies no general singular quotient
classification, asymptotic error bound or interval of nonvanishing denominator.

The optional F1–F4 findings require no new proof scope: fresh paired controls
use canonical witnesses, the physical real-current premises remain explicit,
and no higher-grade agreement follows from F1–F4's grade-zero checks. The user
authorized this fixed P1–P4 packet for Claude/Grok; both reviews are complete.
No commit is authorized.


## Domain coverage and evidence

No positivity, Hermiticity, nonzero baseline or physical channel-completeness
assumption is hidden in the polynomial identities. `Path` permits any index
type; contractions require a finite index type, including the empty type.
The real parameters t, r and epsilon may be zero or negative. Complex amplitude
and matrix coefficients are unrestricted. P3 holds in any commutative field.
There is no new dimension or mode census in this algebraic contract.

| Cases/claim | Formal evidence | Controls and compact native evidence |
|---|---|---|
| Independent eta/sigma and any real path ratio, including zero | `rectangle_decomposition`, `rectangle_path`, `baseline_free`, `path_does_not_identify_rectangle` | `path_ratio` directly checks rectangle evaluation; `delta_zero_jet` checks the delta. The path identity and native wrong mixed-ratio-power control supply the path link. The kernel identity alone does not assert a nonzero kernel vector in an empty amplitude space. |
| Every supplied degree-two amplitude/current coefficient, no sign restrictions | `quadratic_exact`, `retained_with_remainder`, `real_flux_expansion` | `first_interference`, `current_variation`, `second_amplitude`, `baseline_interference`, `higher_terms`; all seven native degrees and full off-diagonal current fixture. |
| j0 nonzero | `quotient_equations`, `quotient_unique`, `quotient_residual` | `leading_normalization`, `incident_first`, `incident_second`, `negative_leading`; complex-field positive and per-column native recurrence/residual checks. |
| j0 zero versus nonzero is exhaustive | `leading_denominator_cases` | The ordinary quotient conclusions explicitly require j0 nonzero. Native nonfinite zero-denominator output is an invalid-domain witness, not valid quotient correspondence. |
| j0=0,n0 nonzero; j0=n0=0 | `zero_leading_obstruction`; `zero_model_two_solutions` supplies distinct constants for the all-zero leading equation | Obstruction and nonuniqueness positives; `zero_model_two_constants` pair. Higher coefficient equations in the general singular case are not classified; no general existence or uniqueness is asserted there. |
| Every real epsilon, with a separately defined fraction domain | `flux_epsilon_squared`, `scaled_denominator_nonzero`, `epsilon_cancels` | `epsilon_squared`; native quadratic scaling. A defined scaled fraction requires nonzero epsilon and nonzero incident denominator. Lean's totalized division identity alone does not enforce that domain. |
| Arbitrary baseline/induced amplitudes; zero baseline specialization | `subtracted_total`, `subtracted_eq_induced_iff`, `baseline_zero_total`, `induced_coefficient`, `induced_low_coefficients` | `induced_leading`, `subtracted_is_not_induced`, `omitted_parent_second`; native total-minus-baseline versus induced fixture. The second-amplitude witness varies retained a2, not an omitted pure parent eta²/sigma² coefficient. The iff states precisely when the real interference difference vanishes. |
| Empty finite amplitude space | General contraction theorems and an explicit empty-amplitude positive | Zero contraction is admissible; no incident nonzero-denominator premise is inferred from it. |

Run7 checked seven fresh objects, 33 selected standard-axiom lists, sixteen
paired mathematical rejections and twenty positive executions (44 records with
native evidence). Three omission pairs share the same q2 positive, so the twenty
executions contain eighteen distinct positive statements. Every accepted false
statement has one `contract_control` diagnostic with one unsolved `False` and
no other error or warning. These are instance statements decided using canonical
witnesses or arithmetic, not source-replacement mutations or independent
rederivations. The compact native run1 has thirty identities/examples, eleven
wrong-formula controls and one invalid-domain witness. Its unchanged hashes were
revalidated in run7, rather than executing it again.

See SCATTERING_BOOKKEEPING_FIDELITY.md and SCATTERING_BOOKKEEPING_VERIFICATION.txt
for provenance, resource limits and the translation boundary. Both independent
source-fidelity reviews are CLEAR; neither reviewer reran Lean or NumPy.
SCATTERING_BOOKKEEPING_FIDELITY_REVIEW.md records optional-note dispositions.
The closure manifest records documentation-only deltas from the unchanged
approved snapshot. Conditional finite-solve sensitivity has not started:
the user's instruction to pause after this review supersedes further work.
