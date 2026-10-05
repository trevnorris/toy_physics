# S10 coverage and remaining work

Updated 2026-09-11. This is a mathematical and implementation coverage map,
not a claim that the supplied physical model has been derived.

## Current scope and Lean completion

[FORMALIZATION_POLICY.md](../FORMALIZATION_POLICY.md) governs this work. The
systematic CAS bridge expansion is stopped; the historical bridge remains
in `archive/pre-cleanup-2026-10-04`, with the compact D3 identity kept in the tree. Outstanding transcript fields are not automatically Lean obligations.

**The scoped S10 Lean contract is complete.** C1–C4 below are closed: the existing
six-family proofs supply the mathematical classification, the compact
action/operator connection is identified, the essential controls pass, and
independent Claude and Grok fidelity reviews both returned CLEAR. See
[the review and disposition record](FIDELITY_REVIEW.md) and
[the build record](CONTRACT_BUILD_VERIFICATION.txt).

The Lean task stops here. The CAS production, comparator/export and paper tasks
below retain separate completion criteria. A future discrepancy that invalidates
the compact fidelity link would reopen that named obligation; it would not
justify resuming systematic transcript translation.

## S10 work contract

**Scope:** the conditional real plane-wave classification of the six supplied
quadratic action families. The contract below assembles existing theorems; it
introduces no new physical premise or requirement to translate more CAS output.
Its status is **complete within the stated Lean scope**. The exact reviewed v1
contract is archived (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_fidelity_contract_v1.md`);
subsequent wording clarifications and status updates are recorded in the review
disposition. No canonical proof or control instrument changed after review.

### Claim and domain

Work on flat whole spacetime with constant real `rho > 0`, `mu > 0`, real
wavevector `k != 0`, real amplitudes and spatial dimension `D >= 2`. The actions
use smooth backgrounds and smooth compact test variations. Plane waves are
`u(t,x) = a cos(k·x - omega t)`. Stationarity is against all admissible compact
variations, not only variations within that ansatz. A relative action avoids
assuming finite total action for a nondecaying plane wave.

For ANISO, the distinguished index `e` is any `Fin D`, `sigma > 0`,
`sigma != 1`, and `sigma` is declared dimensionless. For XCOEF_SCALE, the
coefficient `c > 0`, `c != 1` is declared dimensionless. Here `c != 1` labels a
nontrivial control and is needed only for distinctness from MAIN; the root and
count formulas also hold at `c = 1`. SIGNFLIP is `c = -1`. The signed spectral variable `z` means
physical squared frequency considered as an arbitrary real number; a negative
`z` is not a real-frequency cosine wave. ANISO's internal normalized variable
is `rho * z / mu`; it must not be confused with physical `z`.

Let `K = |k|^2`, `R = (mu/rho) K`, `T = {a : k·a = 0}`, and `L = span{k}`.
For ANISO let `p = k_e`, `q = K - p^2 >= 0`,
`E = (mu/rho) (p^2 + q/sigma)`,
`O = {a : a_e = 0 and k·a = 0}` and
`w = (q + sigma p^2) e - sigma p k` (here `e` in the vector expression denotes
its coordinate unit vector). N2 is full kernel dimension; N3 is the dimension
of its intersection with `T`. Neither is a count of displayed basis vectors.

| Package | Nonzero squared-frequency candidate and full amplitude space | N2 / N3 there | Static amplitude space; N2 / N3 |
|---|---|---|---|
| MAIN | `R`, space `T` | `D-1 / D-1` | `L`; `1 / 0` |
| XFORM_FULLGRAD | `R`, whole amplitude space | `D / D-1` | zero space; `0 / 0` |
| XFORM_DIVONLY | `R`, space `L` | `1 / 0` | `T`; `D-1 / D-1` |
| XFORM_SIGNFLIP | `-R`, space `T`; negative spectral branch | `D-1 / D-1` | `L`; `1 / 0` |
| XCOEF_SCALE | `c R`, space `T` | `D-1 / D-1` | `L`; `1 / 0` |
| XFORM_ANISO | `R` and `E`, with the directional cases below | as below | `L`; `1 / 0` |

Positive branches admit real frequencies; SIGNFLIP has no nonzero stationary
real-frequency cosine wave at nonzero real frequency. For FULLGRAD, zero is
not a determinant root on this domain. Elsewhere, frequencies outside the
listed branches have no nonzero amplitude. Algebraic multiplicity is distinct
from these dimensions; the existing imported D3 polynomial proofs certify
multiplicities on their declared scope, not on every possible emitted package.

| ANISO direction | Predicate | Ordinary branch | Extra branch |
|---|---|---|---|
| Parallel | `q = 0` | one merged root `R=E`, space `O=T`, `D-1 / D-1` | same root, counted once |
| Perpendicular | `p = 0`, `q > 0` | `R`, space `O`, `D-2 / D-2` | `E`, space `span{w}`, `1 / 1` |
| Oblique | `p != 0`, `q > 0` | `R`, space `O`, `D-2 / D-2` | `E`, space `span{w}`, `1 / 0` |

These cases are exhaustive and disjoint on `k != 0`: `perpSq_nonneg` and
`perpSq_zero_iff` identify the parallel case, and the remaining case splits
on `k_e = 0`. `frequency_coincidence_iff` proves exactly `q=0` for the merger;
`extra_exactly_transverse_iff` proves exactly `p=0` for the split extra mode's
transversality. A generic chart assumption is not a hypothesis of this census.
At D=2 the ordinary split candidate has zero nullity and is not an actual mode.
The separate D=1 theorems are additional evidence outside this contract's domain.

### Compact fidelity link

The common action is `rho/2 * v^T W v - mu*c/2 * S(J)`. MAIN uses `W=I`,
`c=1` and `S=(1/2) sum_ij (J_ij-J_ji)^2`; FULLGRAD uses
`S=sum_ij J_ij^2`; DIVONLY uses `S=(sum_i J_ii)^2`. SIGNFLIP sets `c=-1`,
XCOEF_SCALE sets `c` to the supplied scale, and ANISO changes only `W_ee` to
`sigma`. All other inertia entries remain one. In Lean the spatial derivative
`J_ij` is `J i.succ j`; `J 0` is velocity. The plane-wave operator is
`rho*z*W - mu*c*B(k)`, with `B=K I-k k^T` for curl stiffness, `B=K I` for
FULLGRAD and `B=k k^T` for DIVONLY.

The actual definitions are selected by
[Packages.lean](S10Audit/Packages.lean), with
`S10Audit.packageAction_uses_stiffness`; this is an action identity, not an
identification inferred from matching spectra. The underlying definitions are
[baseline Action](S10Pilot/Action.lean), [form-control Action](S10Controls/Action.lean),
[Scalar](S10Controls/Scalar.lean) and [anisotropic Action](S10Anisotropic/Action.lean).

The CAS source connection is the six-entry `PACKAGES` selector,
`stiffness_density` and `build_action` in
[the SymPy builder](../../scripts/S10_brane_mode_spectrum_sympy_audit.py), and
`buildPackage` in [the Wolfram builder](../../mathematica/S10_brane_mode_spectrum_mathematica_audit.wl).
Map `rho_br` / `rhoBr` to `rho`, `mu_R` / `muR` to `mu`, `s_rho` / `sRho` to
`sigma`, and `s` / `coefficientScale` to `c`. Both CAS engines distinguish their
first displacement/velocity component (the amplitude index); this is Lean index
`0`. Source correspondence for all
six actions is reviewed at this compact construction boundary. It is not a
claim that Lean executes or certifies either builder.

The stronger existing ANISO D3 artifact link is retained:
[Bindings.lean](S10Audit/CAS/Bindings.lean) proves each imported route-B matrix
is one half of the same action-derived operator, with `M_A = -2 M_B`.
The scalar factors are nonzero and preserve the kernel, which is the claim
here. This does not establish equality of normalized resolvents or residues.
The focused comparator (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/S10_anisotropic_strata_comparator.py`)
checks both exceptional directions at D=3,4, actual bases and normalized roots;
its result (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/scripts/out/S10_anisotropic_strata_comparator.json`) reports
sampled coverage. General coverage comes from the theorem census above.

All six actual action constructors and the limited stronger D3 link are in the
review packet. Remaining broad-comparator parsing, naming and sign findings do
not become proven by this contract. Any finding that changes the action/operator
identification must be resolved before the fidelity link is cleared.

### Existing theorem and control evidence

All names below are qualified by their file's namespace. No new standalone
proof is required merely to combine existing checked statements in this table.

| Contract obligation | Existing theorem anchors | Mutation / admissibility evidence |
|---|---|---|
| Action, variation and actual PDE; phase-average normalization | MAIN, form controls and ANISO: `lagrangian_variation`, `actionStationary_iff_eulerLagrange`, `eulerLagrange_planeWave`, `phaseAverage_eq`; scalar inherits MAIN stationarity via `lagrangian_eq`, `actionStationary_eq`, with `variational_planeWave_iff` and `phaseAverage_eq` | original action/form mutations in RESULT, CONTROLS_RESULT and ANISOTROPIC_RESULT; compact closure checks reproduce these with passing originals |
| MAIN complete census and S9 specialization | `S10Pilot.s10_variational_certificate`, `propagating_variational_mode_iff`, `S10Pilot.Specialization` equalities; static N3 from `S10ScalarControls.zero_counts` at `c=1` and `S10Controls.longitudinal_inf_transverse` | nonzero wave existence and longitudinal counterexamples; stiffness normalization mutation |
| FULLGRAD and DIVONLY full/static spaces | `S10Controls.stiffness_comparison_certificate`, `full_variational_iff`, `div_propagating_variational_iff`, `full_mode_counts`, `div_mode_counts`, `full_zero_space`, `div_zero_space`; DIVONLY static N3 follows from `T ∩ T = T` | wrong stiffness-form mutations; `longitudinal_control_witness`, `transverse_control_witness` |
| Coefficient and sign controls | `S10ScalarControls.nonzero_root_iff`, `cone_modeSpace`, `zero_modeSpace`, `positive_variational_certificate`, `negative_control_no_real_wave`, `negative_root_exists`, `coefficient_changes_frequency` | ignored coefficient/sign mutations and admissible concrete controls; `coneValue_scaling_ratio` |
| Exhaustive ANISO cases and full subspaces | `S10Anisotropic.split_propagating_iff`, `parallel_kernel_iff`, `zero_variational_iff`, `split_variational_certificate`, `oblique_counts`, `perpendicular_counts`, `parallel_counts`, `zero_counts` | existing rejected `missing_parallel_branch`, `missing_perpendicular_branch`, `incomplete_parallel_kernel`, `dependent_parallel_basis`, `generic_count_beyond_chart`; concrete parallel/perpendicular/oblique controls |
| Root coincidence, positivity and multiplicity distinctions | `frequency_coincidence_iff`, `extra_exactly_transverse_iff`, `nonzero_root_positive`; existing D3 RootBindings on their stated scope | rejected `collapsed_parallel_multiplicity`, `incomplete_cubic_root_multiset`, `confused_candidate_and_distinct_count`, `wrong_extra_root_sign`, `positive_root_without_nonzero_wavevector` |
| Q5/Q6/Q7 mathematical assertions | scalar and ANISO scaling theorems; `S10Audit.Expr.HasDim.rescale`, dimension solves, root dimensions, `packageAction_uses_stiffness`, `epsilonCurl_eq_jetCurl`, six Q7 comparisons | existing rejected root-wavevector, slot-unit and epsilon-contraction mutations; `q7_control_counterexample`, coefficient free-family and declared-dimensionless theorems |
| Actual selected D3 matrix link and denominator domain | `S10Audit.CAS.PY/WL.matrixB_action`, `matrixB_reference`, `routes` | rejected PY/WL wrong coefficients, `missing_chart_denominator`, `wrong_basis_normalization`; passing denominator and basis controls |

The durable existing records are CAS controls (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_cas_bridge_checks.json`),
Q6/Q7 controls (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_q6_q7_checks.json`),
matrix controls (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_matrix_checks.json`) and
CAS verification (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/CAS_BRIDGE_VERIFICATION.txt`). The [closure instrument](../../_measurements/S10_lean_contract_check.py) and
[its result](../../_measurements/S10_lean_contract_checks.json) reproduce the
early action controls and add scalar/sign and phase-average controls: seven
mathematical rejections and six passing controls, with commands, source variants
and diagnostics retained. Canonical proof files are unchanged.

### Review, exclusions and finite completion

The author is Codex. Fresh Claude and Grok reviews independently inspected the
same fixed v1 contract and source/evidence packet; both returned CLEAR with no
substantive blocker. They inspected definitions, premises, quantifiers,
normalization, excluded cases, correspondence and mutation meaning. The
[review record](FIDELITY_REVIEW.md) retains identities, packet hash, reports,
editorial dispositions and limits. Reviewers did not rerun the build or controls.

Completion consists only of: (C1) this contract mapped to existing theorems;
(C2) the compact action/operator fidelity connection checked at the stated level;
(C3) the essential controls documented with passing counterparts and substantive
failure evidence; (C4) both independent reviews resolved, with build/axiom and
source-provenance evidence retained. Only an identified gap in C1-C4 justifies
further Lean implementation. **All four completion items are met.**

Excluded: selecting physical D=3 or deriving the action, arbitrary Fourier/PDE
completeness, interfaces, dissipation, nonlinear/strained backgrounds, an explicit
exponentially growing spacetime construction, general CAS discovery correctness,
and whole-transcript certification. Production Q7, CAS stratum handling, the
broad comparator/export refresh and ledger/paper reconciliation remain separate
S10 tasks. This contract does not declare the full ledger pipeline complete.

## Existing coverage and remaining obligations

| Requirement | Current coverage | Remaining obligation |
|---|---|---|
| Supplied action and actual first variation | Baseline plus all five controls in Lean | Physical justification remains a premise. |
| Compact-test stationarity and local PDE | All six action families | Boundary/interface and weak-solution extensions are outside this S10 setting. |
| Real cosine phase average | All six actions; exact 0..2pi integral; actual anisotropic D3 matrices now prove `M_A = -2 M_B` and their action linkage | Compact normalization link complete above. Production route comparison remains separate; do not replicate the full expression bridge. |
| Baseline and form-control spectrum | Complete real plane-wave amplitude classification in arbitrary D | This does not prove completeness of general PDE solutions. |
| Anisotropic spectrum and exceptional directions | Complete oblique, parallel, perpendicular, and static classifications in Lean | The focused CAS samples do not certify a general stratum-discovery algorithm. |
| Coefficient and sign controls | Arbitrary real squared-frequency classification; positive/negative sign and N2/N3 counts | An exponential spacetime solution is not yet constructed. That construction is outside the current claim and requires an explicit contract extension. |
| Q5 scaling | Nonzero root formulas and ratios, including the anisotropic extra branch | Zero-root ratio remains undefined, rather than assigned a value. |
| Q6 dimensions | PhysLean expression checker and unit-change theorem; actual-action bridges, coefficient solve/inventory, full root ratios, Q7, matrices, mixed-unit minors, complete basis families and N5/N6 residuals; kernel/rank unit invariance and vacuity proved. Both engines' generic and exceptional anisotropic D3 matrices, determinants, roots, stacks, bases, residuals and emitted minors have checked expression bridges (916 scalar expressions including generic/exceptional counts, root lists, aggregate coincidence differences and spectrum solve operands), explicit domains and complete-kernel basis proofs. Fixed-point reruns have an explicit coordinate/physical scale convention. | Necessary dimensional/scaling invariants and parameter map are included in the completed contract. Remaining per-output metadata and other emitted cases belong to CAS/comparator coverage; no systematic Lean expansion is planned. |
| Q7 ordinary three-dimensional curl | PhysLean Levi-Civita contraction and all six package comparisons, with exact action-stiffness linkage | Align the CAS implementations and production comparator with the explicit construction. |
| Q8 focused CAS repair | Both engines inspect N2 and N3 matrices; ten sampled roots compare successfully. Lean now checks all 83 emitted D3 minors, completeness of the row/column selections, all 12 rank-drop locus predicates and all four targeted points, including both branches of the extra transverse locus. The rerun matrices, roots and all 12 printed basis vectors are connected to complete kernel proofs and physical rescaling. All 112 generic and exceptional rank/nullity/basis-count records are certified against the actual objects; generic bindings keep explicit chart assumptions and N4/N7 use signed subtraction. All six solution/distinct-root list pairs have determinant-completeness proofs; distinct counts, algebraic multiplicities and Wolfram candidate-filter counts are certified. The primary coincidence equations, guarded loci, allowed regions, Boolean/outcome decisions and SymPy witnesses are certified on the positive-coefficient domain, including the full parallel axis. All Q8/aggregate coincidence fields, all 16 reported root signs, eight empty root solver-condition lists, three spectrum solve operands/statuses and six retained/skipped-stratum records are now bound to their mathematical meanings. Lean resolves both SymPy undecided extra-root signs as positive. | Exhaustive cases and invariants are consolidated in the completed contract above. Production stratum handling and full comparator/export integration remain separate tasks; further metadata translation into Lean is not required. |
| Ledger/paper alignment | New evidence recorded in S10 and linked proof reports | Reconcile the paper's historical claims and evidence pointers with the completed formal coverage. |

The focused CAS rerun and its comparison are documented in
[S10_anisotropic_strata_report.md](../../_measurements/S10_anisotropic_strata_report.md).
The latest control results are in [SCALAR_RESULT.md](SCALAR_RESULT.md).
The dimensional-analysis and curl extension is in [Q6_Q7_RESULT.md](Q6_Q7_RESULT.md).
The matrix, minor and complete-basis extension is in [MATRIX_RESULT.md](MATRIX_RESULT.md).
The first actual-expression bridge is in [CAS_BRIDGE_RESULT.md](CAS_BRIDGE_RESULT.md),
with current build evidence in CAS_BRIDGE_VERIFICATION.txt (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/s10/CAS_BRIDGE_VERIFICATION.txt`).
The complete minor and exceptional-locus bridge is in [MINOR_LOCUS_RESULT.md](MINOR_LOCUS_RESULT.md).
The exceptional matrices and complete bases are in [EXCEPTIONAL_RERUN_RESULT.md](EXCEPTIONAL_RERUN_RESULT.md).
The exceptional rank/nullity and signed-count bindings are in [COUNT_RESULT.md](COUNT_RESULT.md).
Root-list completeness and multiplicities are in [ROOT_RESULT.md](ROOT_RESULT.md).
The primary coincidence loci, decisions and witnesses are in [COINCIDENCE_RESULT.md](COINCIDENCE_RESULT.md).
The aggregate/Q8 fields, root signs, solve operands and stratum dispositions are in [METADATA_RESULT.md](METADATA_RESULT.md).

S11c-d currently consumes the frozen S11c-b base and S11c-c1/c2 deltas.
This increment preserves those files and `S10_exports.py`; it makes no changes
to the S11c-d constant-end resolvent or scattering construction. Any proposed
change to those shared physical operands requires a separate dependency review.
