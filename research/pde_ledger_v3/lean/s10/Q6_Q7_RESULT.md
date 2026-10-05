# S10 dimensions and Levi-Civita comparison

2026-09-11. The `S10Audit` library extends the existing S10 action proofs with
expression-tree dimensional analysis and the ordinary three-dimensional curl
comparison. It uses the pinned PhysLean dimensional algebra and Levi-Civita
symbol directly. No supplied physical premise is derived by these checks.

## Q6: action expressions, coefficient solve and roots

[Dimensions.lean](S10Audit/Dimensions.lean) defines an expression language for
atoms, real constants, addition, subtraction, multiplication, division, natural
powers and nonempty finite sums. The dimension checker walks every node and
requires matching dimensions on every additive branch. Its successful output
is equivalent to a dimensional typing derivation. A second theorem proves that
every such derivation transforms correctly under multiplicative changes of
units. These are properties of the entire expression language, not only the
particular S10 formulas.

[ActionTrees.lean](S10Audit/ActionTrees.lean) represents each supplied action,
including every field-derivative factor and coefficient. Evaluating each tree
is proved equal to the existing action definition for arbitrary real inputs
and arbitrary jets. Bare displacement atoms are also supported and retain
their field dimension; a squared bare displacement cannot silently acquire
dimensionless units. The action tree bridge checks algebraic equality to the
Lean density; it does not parse or certify a CAS transcript.

The solve starts with an arbitrary displacement dimension `U`. In multiplicative
dimension notation, the inferred coefficient dimensions are

```text
rho = energyDensity / (U / time)^2
mu  = energyDensity / (U / length)^2.
```

[DimensionSolve.lean](S10Audit/DimensionSolve.lean) derives these relations from
the inferred dimensions of the actual kinetic and stiffness action-term trees.
For the ordinary-inertia packages, these equations characterize all solutions
for the coefficients that occur in those trees. The physical coefficient census
here uses D>=2, matching the S10 sweep. The checker retains cancelling terms
and checks structural homogeneity; it does not minimize constraints by algebraic
cancellation, such as the identically zero curl stiffness at D=1.
Without the scale declaration, the coefficient
control determines only `scale * mu`; an explicit theorem retains the complete
free family obtained by choosing any scale dimension. The supplied
`scale = dimensionless` condition then fixes `mu`.

For one-axis anisotropic inertia in D>=2, an unweighted kinetic component fixes
`rho`, and the distinguished kinetic component then forces `sigma` to be
dimensionless. No dimensionless premise for `sigma` enters that theorem. The
D>=2 condition matters: a second inertia component is needed for this argument.

Specializing the displacement premise to length gives the exact L/T/M vectors
for arbitrary D:

| Coefficient | Derived dimension |
|---|---|
| `rho` | `(-D, 0, 1)` |
| `mu` | `(2-D, -2, 1)` |
| anisotropic `sigma` | `(0, 0, 0)`, derived for D>=2 |
| coefficient `scale` | `(0, 0, 0)`, supplied declaration |

[CoefficientInventory.lean](S10Audit/CoefficientInventory.lean) traverses each
action tree to identify its coefficient symbols. The declared scale is removed
from the unknown set; the solved anisotropic scale stays in it. This gives two
unknown dimension vectors for the ordinary-inertia packages and three for
anisotropy, or six and nine scalar L/T/M exponents respectively. The solve
theorems characterize the corresponding independent constraints; the
anisotropic three-coefficient block is checked to have determinant one.

The action homogeneity result is explicitly a consequence of the solve for
**every** supplied `U`. It therefore cannot independently verify `U`. A separate
theorem proves that `mu/rho` has dimension `(length/time)^2` regardless of `U`.
Another proves that multiplying an action by any dimensionless numerical
constant leaves the checker green. These limitations are retained as results.

[RootDimensions.lean](S10Audit/RootDimensions.lean) constructs the entire
ordinary, coefficient/sign and anisotropic-extra root expressions. Their real
evaluations equal the already verified cone formulas. For the solved units,
the trees have dimension `time^-2`, and dividing them by the full wavevector
norm-squared tree gives `length^2/time^2`. The extra branch includes its
direction-dependent numerator and inertia divisor. A negative squared-frequency
branch has the same units as a positive one; dimensions do not test its sign.

Dropping the wavevector factor is proved to change the inferred root dimension.
No unique physical unit is assigned to the identically zero root here. The
literal-zero syntax has the dimensionless convention used for numerical
constants; algebraically equal zero expressions can carry other dimensions.

## Q7: the exact stiffness operand and ordinary curl

[Packages.lean](S10Audit/Packages.lean) selects the existing action and the exact
stiffness definition it uses, with a theorem connecting the two. The stiffness
coefficient and its sign remain separate from that operand, as required by Q7.

[Curl.lean](S10Audit/Curl.lean) defines the ordinary curl by
`c_i = sum(j,k) epsilon_ijk J[j+1,k]`, using PhysLean's determinant-based
Levi-Civita symbol. Only after that construction is it proved equal to the S9
component definition. The resulting squared norm equals the baseline
antisymmetric double-sum stiffness with its supplied one-half normalization.

| Package stiffness | Proved difference from ordinary curl-squared |
|---|---|
| MAIN, SIGNFLIP, ANISO, XCOEF_SCALE | zero |
| FULLGRAD | `sum(i,j) g_ij g_ji` |
| DIVONLY | `(sum(i) g_ii)^2 - sum(i,j) g_ij^2 + sum(i,j) g_ij g_ji` |

A symmetric diagonal gradient gives zero curl-squared and difference one for
each form control. The Q7 expressions in [Checks.lean](S10Audit/Checks.lean)
also evaluate to these exact operands and receive a proved whole-expression
dimension. Q7 remains a form/normalization comparison, not a mode-count proof
or a test of the sign with which stiffness enters the action.

## Verification and remaining boundary

Run from `research/pde_ledger_v3/lean`:

```sh
LAKE_CACHE_DIR=.lake/cache lake build
```

The committed checkpoint build succeeded with **250 audited declarations**, including 99 in
`S10Audit`, and only `propext`, `Classical.choice` and `Quot.sound`.
There are no proof admissions. Two isolated source mutations are rejected:
doubling the epsilon contraction breaks its curl identity, and dropping the
wavevector factor breaks both the root-evaluation and root-dimension theorems.
The mutation instrument (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_q6_q7_check.py`) and
recorded checks (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/_measurements/S10_lean_q6_q7_checks.json`) retain that evidence.

The local cache setting permits newly used PhysLean modules to build in the
workspace when the toolchain installation is read-only. Dependency pins are
unchanged. [AUDIT_VERIFICATION.txt](AUDIT_VERIFICATION.txt) records the theorem
audits, build result, source hashes and targeted mutation checks.

This completes the Lean action/root core of Q6 and the six-package Lean Q7
comparison. The later [matrix and basis extension](MATRIX_RESULT.md) adds
expression-tree dimensions for action-derived matrices, mixed-unit stacked
minors, complete coordinate bases and N5/N6 residuals, with kernel/rank
invariance under unit changes. The subsequent [CAS bridge pilot](CAS_BRIDGE_RESULT.md)
connects the actual generic anisotropic D3 expressions from both engines. The
[exceptional extension](EXCEPTIONAL_RERUN_RESULT.md) adds the targeted rerun
matrices, roots, full bases and explicit coordinate-to-physical scale proofs.
The [count bridge](COUNT_RESULT.md) further certifies the generic and exceptional
N2/N3 ranks and nullities and the signed N4/N7 residuals, with explicit chart
assumptions on the generic records. The [root-list extension](ROOT_RESULT.md)
adds complete solution lists, distinct-root counts, algebraic multiplicities
and syntactic filter counts. The [coincidence extension](COINCIDENCE_RESULT.md)
adds primary root differences, guarded loci, allowed regions, decisions and witnesses.
The [metadata extension](METADATA_RESULT.md) connects Q8/aggregate fields, root
signs, solver operands/statuses and stratum dispositions.
Other emitted cases and the production CAS/export refresh remain open.
The original CAS Q7 implementations still
require alignment with their explicit Levi-Civita construction requirement.
The [coverage map](COVERAGE.md) keeps those tasks separate. These results are
included in the combined S9/S10 checkpoint (`archive/pre-cleanup-2026-10-04:research/pde_ledger_v3/lean/CHECKPOINT.md`).
