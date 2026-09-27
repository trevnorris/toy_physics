# Conditional finite-solve sensitivity: S1–S4

Authorized by “Great. Let's continue.” after the finite-solve sensitivity
recommendation, 2026-09-26. Governed by ../FORMALIZATION_POLICY.md.
Status: local verification PASS in recorded run6: eight fresh isolated objects,
51 selected standard-axiom audits, eighteen paired mathematical rejections and
twenty-four positive executions (twenty-three distinct statements). The compact
native record passed 36 checks in run1 and was hash-revalidated in run6.
Claude and Grok independently returned CLEAR with no blocking findings.
Bounded S1–S4 is complete; optional dispositions and documentation-only closure
are recorded in SCATTERING_SENSITIVITY_FIDELITY_REVIEW.md. F1–F4 and P1–P4
remain complete, with their evidence unchanged.

## Supplied finite object

The actual native anchor is `_measurements/S11c_d_finite_scattering.py:construct`.
Its square collocation matrix M includes the replaced modal-boundary rows; b
contains the prescribed incident data. Row scaling R and unknown scaling D give
K = R M D⁻¹, g = R b, z = D c. The approximate solve uses least squares and then
unscales c. The recorded scaled residual is R(M c - b) = K z - g. A measured
balanced condition number is not a certified bound on ‖K⁻¹‖, and the residual
of the unreplaced bulk matrix is not the residual of this boundary-value system.

For a fixed incident column, trace evaluation T, outgoing trace basis V and
open-channel selection P give a = P V⁻¹(T D⁻¹ z - t_in). Thus observation is
affine, a = C z + d, with the prescribed incoming trace subtracted. The fixed
offset cancels in a difference, but not in the amplitude about which the
quadratic current is evaluated. Native J is restricted to selected open outgoing
channels, with same-end blocks and structurally zero cross-end entries.
The complete supplied outgoing matrix J gives
q(a) = Re(conj(a)ᵀ J a). The native numerator retains off-diagonal entries;
the incident denominator is the signed incoming current, separately required
positive in the native constructor. Origin rephasing transforms both amplitude
and metric, as in the already reviewed F1–F4 pullback identity.

Norms and field/equation unit frames are premises. The operator theorem is
stated on complex normed spaces with a supplied continuous linear equivalence;
finite-dimensional realizations specialize it. It does not certify invertibility
from floating-point rank or singular values. Amplitude bounds use the explicit
finite ℓ¹ mass ∑ᵢ|aᵢ|; a component observation map Cᵢ has bound ∑ᵢ‖Cᵢ‖.
The current coefficient bound β is an entrywise absolute bound, giving
|conj(x)ᵀ J y| ≤ β mass(x) mass(y). This is not silently identified with the
native Euclidean condition number or a spectral matrix norm.

## Bounded obligations

| Item | Claim and cases | Controls |
|---|---|---|
| S1 | Exact residual identity and ‖z_hat - K⁻¹g‖ ≤ κ‖K z_hat - g‖ for a supplied inverse bound. Show row/unknown scaling correspondence. Reuse T2 inverse stability for a supplied perturbation E with κε<1, adding numerical residual and right-hand-side error explicitly. | Wrong residual sign, omitted κ, wrong scaling/unscaling, false small-residual implication; exact-solve and zero-data positives. |
| S2 | Affine observation error bound with fixed C,d, including empty amplitude spaces. Bound the full F1–F4 current change by β(2 mass(a) mass(e) + mass(e)²), and separately account for a supplied current-matrix change. Compose with residual and observation estimates. | Omitted observation factor, baseline offset, interference term, quadratic error or current change; complex/off-diagonal fixtures and zero error positives. |
| S3 | Propagate real numerator and denominator errors into a normalized fraction under explicit positive absolute-denominator lower bounds; prove a denominator margin from a reference value and smaller error. Allow signed numerators/denominators, but do not infer positivity or conservation. State the zero-denominator exclusion and small/zero-margin limitations. | Omitted denominator error or lower-bound factor, unsafe zero denominator; positive and negative admissible denominator examples. |
| S4 | Compact source/operation identification, standard-axiom audit, true/false controls and two independent fidelity reviews of a fixed packet. All claims conditional on the supplied finite operator and maps. | Canonical witnesses/arithmetic plus small selected native AST operations; compiler/resource/import errors never count as rejected mathematics. |

The proof/control map is explicit:

| Claim | Canonical evidence | Paired control names / admissible cases |
|---|---|---|
| S1 residual, scaling and perturbation | Residual.residual_identity/residual_bound, exact_residual/residual_zero_iff, balanced_apply/balanced_residual, unscale_error, perturbed_residual_bound; reviewed T2 inverse_norm_bound/solution_error_bound | residual_sign, inverse_factor, small_residual_large_error, row_column_scaling, unscale, critical_inverse_margin; actual_residual_bound_positive and actual_zero_residual_positive |
| S2 full current and affine observation | Current.pair_bound, flux_bound, flux_change_bound/flux_change_budget, current_change_bound/flux_and_current_change, amplitude_difference/amplitude_error_bound | observation_factor, affine_baseline, cross_term, quadratic_error, current_change, full_complex_current; empty_amplitude_positive and actual_flux_bound_positive |
| S3 normalized fractions and margins | Fraction.denominator_margin/denominator_nonzero, fraction_difference, fraction_error_bound/fraction_error_budget, denominator_cases | normalized_error, denominator_change, denominator_lower_bound, negative_denominator, denominator_margin, zero_margin; actual_denominator_margin_positive and actual_signed_fraction_positive |
| Composition | Pipeline.observed_residual_bound/observed_perturbed_residual_bound, flux_distance_budget/flux_residual_bound, normalized_residual_bound | The preceding scale, observation, current and fraction controls target each factor; selected native synthetic conditional-bound examples test their interpretation. |

Module prefixes above identify files; declarations share the namespace
S11ScatteringSensitivity. Paired controls are fresh instance-level statements
resolved with canonical witnesses/arithmetic. They are not mutations of the
canonical proof source or independent derivations. The two current-error
omission controls share one positive statement. The critical-margin pair tests
the zero scalar coefficient; critical_inverse_margin separately proves that
this scalar equation has no inverse solution. No inverse-existence conclusion
at the critical boundary is inferred from the arithmetic pair alone.

Coverage uses separate axes, not a count of physical strata. Exact/nonexact
residual cases are characterized by residual_zero_iff. Empty/nonempty finite
amplitude spaces are both admitted; the empty case has mass and current zero.
Fraction denominators split disjointly into zero and positive absolute value
(denominator_cases), covering either sign. Existing abs_pos/abs_zero give
disjointness; no separate local named disjointness lemma is claimed. Estimates
require positive margins;
a zero margin can admit a zero denominator. The perturbation estimate is
conditional on κε<1; κε≥1 is an excluded estimate domain, not a theorem that
every such operator is singular. The critical scalar witness shows why a
universal extension would fail. These axes may overlap independently.

All real error budgets are nonnegative where needed. Finite amplitude index
types may be empty. No Hermitian, positive-current or nonzero-baseline premise
is imposed on the quadratic algebra. The inverse premise excludes a singular
finite solve; no inverse estimate is asserted at κε≥1. The fraction is interpreted
only on its explicit nonzero domain, irrespective of Lean's totalized division.

## Fidelity, verification and stopping rule

Read-only anchors also include `S11c_d_finite_scattering_resolution.py:inspect`,
`S11c_d_finite_scattering_domain.py:phases/observable`, the shared physics and
scattering-form amendment. Selected original arithmetic/NumPy AST statements
may execute only on small synthetic matrices and traces, with runtime versions,
source/AST hashes and numerical tolerances recorded. No production module
import, saved scientific operand, physical solve, S11c output or driver rerun.
Native source translation stays outside the Lean kernel. Existing empirical
finite-resolution comparisons remain empirical, not newly certified budgets.

New isolated proof objects; every historical source/report and old object is
preserved. Scientific jobs use scripts/s11c_guarded_run.py, one at a time,
2 GiB/no swap/one CPU/32 tasks, Lean -j1 -M4096 with strict warnings and
180-second per-process limits. The internal Lean accounting setting was raised
after run1 hit its 2048 MiB interpreter limit while loading unchanged Stability;
the whole-job cap remains 2 GiB. Actual run1 peak was 1,146,425,344 bytes with
no max/OOM/swap events. No unguarded fallback or automatic retry.
Use silent local completion/error hooks without healthy-job polling.

Local proof, correspondence, controls, audit and preservation obligations passed
at the recorded scope, and both independent reviews cleared. The user authorized
this fixed packet transfer; the approved snapshot remains unchanged. The closure
record documents only live wording/status changes. No commit is authorized.
Portable-runner registration is separate tooling work, not a claim of this contract.

Excluded: constructing physical inverse bounds, certified numerical intervals,
full-operator discretization/quadrature/tail errors, transparent boundaries,
continuum scattering convergence, physical channel completeness, parent-theory
Taylor accuracy, conservation, new pole searches or a systematic CAS bridge.
