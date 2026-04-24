#!/usr/bin/env python3
"""
Stage V2-12: STF angular source-map theorem audit.

This script verifies, symbolically, that the canonical real STF l=2
angular basis gives an identity angular source map between an orbital/worldtube
STF quadrupole and the grouped real P2 throat ports.

The point of the audit is to separate the angular question from the still-open
radial/axial and outgoing-normalization questions.
"""

import sympy as sp

sqrt = sp.sqrt
pi = sp.pi

def tr_inner(A: sp.Matrix, B: sp.Matrix) -> sp.Expr:
    """Euclidean STF tensor inner product Tr(A B)."""
    return sp.simplify((A * B).trace())

def matrix_zero(M: sp.Matrix) -> bool:
    return all(sp.simplify(x) == 0 for x in list(M))

def vector_zero(v) -> bool:
    return all(sp.simplify(x) == 0 for x in list(v))

# ---------------------------------------------------------------------------
# 1. Canonical real STF basis.
# ---------------------------------------------------------------------------
E20 = sp.diag(-1/sqrt(6), -1/sqrt(6), 2/sqrt(6))
E21c = sp.Matrix([[0, 0, 1/sqrt(2)],
                  [0, 0, 0],
                  [1/sqrt(2), 0, 0]])
E21s = sp.Matrix([[0, 0, 0],
                  [0, 0, 1/sqrt(2)],
                  [0, 1/sqrt(2), 0]])
E22c = sp.diag(1/sqrt(2), -1/sqrt(2), 0)
E22s = sp.Matrix([[0, 1/sqrt(2), 0],
                  [1/sqrt(2), 0, 0],
                  [0, 0, 0]])

labels = ["20", "21c", "21s", "22c", "22s"]
basis = [E20, E21c, E21s, E22c, E22s]

trace_checks = [sp.simplify(E.trace()) for E in basis]
gram_stf = sp.Matrix([[tr_inner(A, B) for B in basis] for A in basis])
gram_stf_expected = sp.eye(5)

# ---------------------------------------------------------------------------
# 2. Harmonic normalization using the exact fourth-moment identity.
#
# Y_A(n) = sqrt(15/(8*pi)) E_A^{ij} n_i n_j.
#
# Since E_A is trace-free,
# int (E_A nn)(E_B nn) dOmega = (8*pi/15) Tr(E_A E_B).
# Therefore int Y_A Y_B dOmega = Tr(E_A E_B).
# ---------------------------------------------------------------------------
harmonic_prefactor = sqrt(sp.Rational(15, 1)/(8*pi))
angular_gram_from_moments = sp.Matrix([
    [sp.simplify(harmonic_prefactor**2 * (8*pi/15) * tr_inner(A, B))
     for B in basis]
    for A in basis
])

# ---------------------------------------------------------------------------
# 3. Optional direct integration check using explicit theta/phi harmonics.
# ---------------------------------------------------------------------------
theta, phi = sp.symbols("theta phi", real=True)
nvec = sp.Matrix([
    sp.sin(theta)*sp.cos(phi),
    sp.sin(theta)*sp.sin(phi),
    sp.cos(theta),
])
Y = [sp.simplify(harmonic_prefactor * (nvec.T * E * nvec)[0]) for E in basis]

# Integrate all entries directly. This is not required for the theorem, but it
# catches normalization slips in the explicit real-harmonic realization.
direct_gram_entries = []
for Ya in Y:
    row = []
    for Yb in Y:
        integrand = sp.trigsimp(Ya * Yb * sp.sin(theta))
        val = sp.integrate(sp.integrate(integrand, (phi, 0, 2*pi)), (theta, 0, pi))
        row.append(sp.simplify(val))
    direct_gram_entries.append(row)
direct_angular_gram = sp.Matrix(direct_gram_entries)

# ---------------------------------------------------------------------------
# 4. Source-map identity.
#
# If S(n) = sum_B S_B Y_B(n), then projected port coefficients are
# P_A = int Y_A S dOmega = S_A.
# ---------------------------------------------------------------------------
S_symbols = sp.symbols("S20 S21c S21s S22c S22s")
Svec = sp.Matrix(S_symbols)
projected_source = angular_gram_from_moments * Svec
source_identity_residual = sp.simplify(projected_source - Svec)

# ---------------------------------------------------------------------------
# 5. Orbital/worldtube STF quadrupole extraction and reconstruction.
# ---------------------------------------------------------------------------
x, y, z, mu = sp.symbols("x y z mu", real=True)
r2 = x**2 + y**2 + z**2
X = sp.Matrix([x, y, z])
Q = sp.simplify(mu * (X*X.T - sp.eye(3)*r2/3))  # STF tensor mu x_<i x_j>
Q_coeffs = sp.Matrix([sp.simplify(tr_inner(E, Q)) for E in basis])
Q_rec = sp.zeros(3)
for coeff, E in zip(Q_coeffs, basis):
    Q_rec += coeff * E
Q_rec_residual = sp.simplify(Q_rec - Q)

# Angular function identity:
# sqrt(15/(8*pi)) Q_ij n_i n_j = sum_A Q_A Y_A.
S_from_tensor = sp.simplify(harmonic_prefactor * (nvec.T * Q * nvec)[0])
S_from_coeffs = sp.simplify(sum(Q_coeffs[i] * Y[i] for i in range(5)))
angular_function_residual = sp.trigsimp(sp.simplify(S_from_tensor - S_from_coeffs))

# ---------------------------------------------------------------------------
# 6. Norm identity and grouped metric.
# ---------------------------------------------------------------------------
source_norm = sp.simplify((Svec.T * angular_gram_from_moments * Svec)[0])
source_norm_expected = sp.simplify(sum(s**2 for s in S_symbols))
source_norm_residual = sp.simplify(source_norm - source_norm_expected)

x20, x21, x22 = sp.symbols("x20 x21 x22", real=True)
full_grouped_vector = sp.Matrix([x20, x21, x21, x22, x22])
grouped_metric = sp.diag(1, 2, 2)
grouped_vector = sp.Matrix([x20, x21, x22])
full_norm = sp.simplify((full_grouped_vector.T * full_grouped_vector)[0])
grouped_norm = sp.simplify((grouped_vector.T * grouped_metric * grouped_vector)[0])
grouped_norm_residual = sp.simplify(full_norm - grouped_norm)

# ---------------------------------------------------------------------------
# 7. Conversion between an older Pi convention and the normalized q convention.
# q20 = sqrt(2/3) Pi20, q21c = Pi21c, q21s = Pi21s,
# q22c = 2 Pi22c, q22s = 2 Pi22s.
# This map has no angular kernel.
# ---------------------------------------------------------------------------
B_Pi_to_q = sp.diag(sqrt(sp.Rational(2, 3)), 1, 1, 2, 2)
det_B = sp.simplify(B_Pi_to_q.det())
rank_B = B_Pi_to_q.rank()

# ---------------------------------------------------------------------------
# 8. Results.
# ---------------------------------------------------------------------------
checks = {
    "all_STF_traces_zero": all(sp.simplify(t) == 0 for t in trace_checks),
    "STF_basis_orthonormal": matrix_zero(gram_stf - gram_stf_expected),
    "angular_gram_from_moments_identity": matrix_zero(angular_gram_from_moments - sp.eye(5)),
    "direct_angular_gram_identity": matrix_zero(direct_angular_gram - sp.eye(5)),
    "source_projection_identity": vector_zero(source_identity_residual),
    "orbital_STF_reconstructs_exactly": matrix_zero(Q_rec_residual),
    "angular_function_reconstructs_exactly": sp.simplify(angular_function_residual) == 0,
    "source_norm_identity": sp.simplify(source_norm_residual) == 0,
    "grouped_metric_matches_full_pair_norm": sp.simplify(grouped_norm_residual) == 0,
    "Pi_to_q_map_invertible": rank_B == 5 and det_B != 0,
}

print("Stage V2-12 STF angular source-map audit")
print("=" * 72)
print()
print("Canonical real STF labels:", labels)
print()
print("STF trace checks:")
for lab, val in zip(labels, trace_checks):
    print(f"  Tr(E_{lab}) = {sp.sstr(val)}")
print()
print("STF tensor Gram matrix Tr(E_A E_B):")
print(gram_stf)
print()
print("Angular Gram matrix from fourth-moment identity:")
print(angular_gram_from_moments)
print()
print("Direct theta/phi angular Gram matrix:")
print(direct_angular_gram)
print()
print("Source-map residual int Y_A sum_B S_B Y_B - S_A:")
print(source_identity_residual)
print()
print("Orbital/worldtube STF quadrupole coefficients Q_A = Tr(E_A Q):")
for lab, coeff in zip(labels, Q_coeffs):
    print(f"  Q_{lab} = {sp.sstr(sp.factor(coeff))}")
print()
print("Q reconstruction residual Q_rec - Q:")
print(Q_rec_residual)
print()
print("Angular function reconstruction residual:")
print(sp.sstr(angular_function_residual))
print()
print("Source norm:")
print("  int S(n)^2 dOmega =", sp.sstr(source_norm))
print("  expected          =", sp.sstr(source_norm_expected))
print("  residual          =", sp.sstr(source_norm_residual))
print()
print("Grouped metric check:")
print("  full five-mode norm with pair representatives =", sp.sstr(full_norm))
print("  grouped metric norm x^T diag(1,2,2) x        =", sp.sstr(grouped_norm))
print("  residual                                      =", sp.sstr(grouped_norm_residual))
print()
print("Pi-to-q convention map:")
print("  B = diag(sqrt(2/3), 1, 1, 2, 2)")
print("  det(B) =", sp.sstr(det_B))
print("  rank(B) =", rank_B)
print()
print("Checks:")
for name, status in checks.items():
    print(f"  {name}: {'PASS' if status else 'FAIL'}")

print()
print("Summary verdict:")
if all(checks.values()):
    print("  PASS: The STF angular source-map has no angular normalization defect.")
    print("  PASS: mhat_ang = 1 in the canonical normalized real-STF basis.")
    print("  PASS: Remaining quadrupole normalization is radial/axial/port-amplitude data.")
else:
    print("  FAIL: At least one angular source-map check failed.")
