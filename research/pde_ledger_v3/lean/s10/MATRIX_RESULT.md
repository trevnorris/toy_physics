# S10 matrix, minor and basis dimensions

2026-09-11. This increment extends the [Q6 action/root proofs](Q6_Q7_RESULT.md)
to the matrix and vector objects used by the mode audit. It adds 55 audited
declarations to `S10Audit`. The supplied action and physical dimensions remain
premises.

## Action-derived matrices

[MatrixTrees.lean](S10Audit/MatrixTrees.lean) constructs the modal matrix for
all six packages. Its entries evaluate exactly to the existing modal operator
applied to coordinate unit vectors. Multiplication by the resulting matrix is
proved equal to that operator on **every** amplitude vector.

The squared-frequency expression is a parameter. Any tree with dimension
`time^-2` gives the same matrix-entry dimension, so the ordinary, sign-flipped,
coefficient-scaled, anisotropic-extra and static root trees are covered. The
static tree uses `0 * frequency^2`, retaining the declared frequency slot's
units. This does not assign a unique physical dimension to the number zero.
The operator-evaluation bridge is stated at `z = omega^2`; the dimensional
theorem itself does not require a nonnegative value of `z`.

For displacement dimension `U`, each matrix entry has dimension

```text
m = energyDensity / U^2.
```

At `U = length`, this is the L/T/M vector `(-D, -2, 1)`.
The appended wavevector row in the N3 stack has dimension `length^-1`, which
generally differs from `m`. Zero matrix entries retain their slot's units by
construction.

This matrix is the coefficient matrix of the already verified modal operator.
It does not equate the two raw CAS routes: the Hessian of the period-averaged
action carries the existing dimensionless factor of one-half. Aligning that
normalization with the emitted CAS records remains an integration obligation.

## Determinants, minors and unit changes

[MinorTrees.lean](S10Audit/MinorTrees.lean) constructs determinants recursively
by cofactor expansion, including the empty determinant, and proves evaluation
equal to Mathlib's determinant. For every explicit selection of rows and
columns, the minor's dimension is the product of the selected row dimensions.
Consequently an order-q modal-matrix minor has dimension `m^q`. An order-q
stack minor containing the wavevector row once has dimension
`m^(q-1) * length^-1`. The formal theorem retains the full row selection and
also handles repeated selections.

[MatrixCovariance.lean](S10Audit/MatrixCovariance.lean) proves that every
multiplicative change of units rescales each row by a nonzero factor. It
preserves the whole matrix kernel, matrix rank and every minor's zero/nonzero
status. Different row factors are permitted, so these results apply to the
N3 stack as well as the modal matrix. No rank-constancy assumption or numerical
witness enters these unit-change results.

## Complete bases with explicit domains

[BasisTrees.lean](S10Audit/BasisTrees.lean) supplies transverse coordinate charts:

```text
b_j = unit_j - (k_j / k_p) unit_p,    k_p != 0, j != p.
```

The chart vectors are linearly independent. A reconstruction theorem writes
**every** transverse vector as their linear combination. A nonzero wavevector
always admits a pivot, but no fixed coordinate is assumed globally nonzero.
For the anisotropic ordinary space, the free coordinates exclude the
distinguished axis as well as the pivot. Each chart vector then belongs to
that space, and every vector in the space is reconstructed. Nonzero
`perpSq` guarantees a pivot away from the distinguished axis.

[NormalizedTrees.lean](S10Audit/NormalizedTrees.lean) supplies the full coordinate
basis and normalized one-dimensional families:

- `k / k_p`, on `k_p != 0`, spans the same longitudinal space as `k`.
- `extraVector / perpSq`, on `perpSq != 0`, spans the same extra-mode space.

Both normalized families are proved linearly independent on those domains.
All these coordinate expressions have dimensionless component trees, with
exact evaluation bridges. The extra normalization is not used on the parallel
locus `perpSq = 0`; the already classified parallel eigenspace uses the full
transverse chart. The existing spectrum theorems identify which of these
spaces belongs to each action and root. This increment supplies explicit
complete vector families for those spaces; it does not certify the basis
normalization chosen by either CAS engine.

## Residuals

[ResidualTrees.lean](S10Audit/ResidualTrees.lean) proves evaluation bridges and
dimensions for arbitrary uniformly dimensioned vector trees. Writing `B` for
their component dimension:

| Object | Dimension |
|---|---|
| Matrix times basis vector | `m * B` |
| N5: matrix times wavevector | `m * length^-1` |
| N6: wavevector dot basis vector | `length^-1 * B` |
| N6: `normSq(k) b - dot(k,b) k` | `length^-2 * B` |

The normalized bases above specialize these formulas to `B = 1`. These proofs
check every component of the full residual expressions. A dimension result
does not imply that a residual vanishes; that is a separate algebraic claim.

## Verification and next boundary

`LAKE_CACHE_DIR=.lake/cache lake build` passes with 250 audited declarations:
21 S9, 34 S10 baseline, 48 controls, 48 anisotropic and 99 in `S10Audit`.
Only the standard logical axioms `propext`, `Classical.choice` and `Quot.sound`
occur. Warnings are errors, and no proof admissions are used.

The [mutation instrument](../../_measurements/S10_lean_matrix_check.py) tests
three isolated copies: an incorrect appended-row unit, a dropped cofactor
sign and an incorrect transverse-basis sign. Results are recorded in
[S10_lean_matrix_checks.json](../../_measurements/S10_lean_matrix_checks.json).
[AUDIT_VERIFICATION.txt](AUDIT_VERIFICATION.txt) records the build, audits and
source hashes.

The next step is the CAS expression/emission bridge: map actual emitted
matrices, minors, bases and residuals into the checked expression language,
retain denominator domains and root substitutions, and reconcile route
normalizations. The production Q7 construction, broad comparator/export
refresh and paper alignment also remain open in [COVERAGE.md](COVERAGE.md).
This matrix increment changes no CAS engine or frozen export. It is included
in the [combined S9/S10 checkpoint](../CHECKPOINT.md).
