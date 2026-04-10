#!/usr/bin/env python3
"""
3pn_grouped_p2_middle_block_test_audit.py

Audit the next exact 3PN grouped-P2 step after the target pack.

Main idea
---------
The grouped-P2 target pack froze the exact 3PN COM data vector and the first
time-local O(omega^2) kinematic scaffold built from the grouped source norms.
This script asks the next sharp question:

    If one applies the *minimal local demotion* needed to make that front-end
    live at formal 3PN order, what COM slot pattern does it generate, and does
    it match the exact GR residual?

The answer is no.  The demoted grouped-P2 scaffold lands in a 3-parameter,
rank-3 subspace of the 9 interaction slots l6..l14 and obeys six exact slot
relations that the solved GR 3PN target violates.
"""

from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def subbanner(title: str) -> None:
    line = "-" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


banner("PART I — EXACT FRONT-END ORDER AND THE UNIQUE LOCAL DEMOTION")

# Symbols and natural virial bookkeeping.
v2, d, U = sp.symbols("v2 d U", real=True)
c20, c21, c22 = sp.symbols("c20 c21 c22", real=True)

# u^2 = v^2 - d^2.
u2 = v2 - d**2

# Direct grouped front-end from the target pack, with 1/r replaced by U.
# The undemoted front-end is of the form U^2 * [degree-3 polynomial in (U,v^2,d^2)].
L_front = sp.expand(
    c20 * d**2 * (3 * u2 - U) ** 2 / 3
    + c21 * u2 * (u2 - d**2 - U) ** 2
    + c22 * d**2 * u2**2
)
L_front = sp.expand(U**2 * sp.expand(L_front))

# Virial weights: wt(U)=1, wt(v^2)=1, wt(d^2)=1.
def virial_weight_monomial(mon: tuple[int, int, int]) -> int:
    # mon = (power of v2, power of d, power of U)
    p_v2, p_d, p_U = mon
    return p_v2 + (p_d // 2) + p_U

poly_front = sp.Poly(L_front, v2, d, U)
weights_front = sorted({virial_weight_monomial(mon) for mon in poly_front.monoms()})
print("Undemoted front-end virial weights =", weights_front)
if weights_front != [5]:
    raise AssertionError("Unexpected undemoted virial weight pattern.")

# Minimal local isotropic demotion by one inverse orbital weight: multiply by r ~ 1/U.
# Up to an overall normalization, the unique local monomial demotion is therefore 1/U.
L_dem = sp.expand(L_front / U)
poly_dem = sp.Poly(L_dem, v2, d, U)
weights_dem = sorted({virial_weight_monomial(mon) for mon in poly_dem.monoms()})
print("Demoted front-end virial weights =", weights_dem)
if weights_dem != [4]:
    raise AssertionError("Unexpected demoted virial weight pattern.")


banner("PART II — EXACT COM SLOT MAP OF THE DEMOTED GROUPED-P2 FRONT-END")

# Standard reduced COM 3PN interaction basis (slots l6..l14).
basis = {
    6: U * v2**3,
    7: U * v2**2 * d**2,
    8: U * v2 * d**4,
    9: U * d**6,
    10: U**2 * v2**2,
    11: U**2 * v2 * d**2,
    12: U**2 * d**4,
    13: U**3 * v2,
    14: U**3 * d**2,
    15: U**4,
}

coeffs = {}
for i in range(6, 15):
    coeffs[i] = sp.Poly(L_dem, v2, d, U).coeff_monomial(sp.Poly(basis[i], v2, d, U).monoms()[0])
    print(f"l{i} =", coeffs[i])

# Slots outside l6..l14 should vanish.
expect_zero("slot l15", sp.Poly(L_dem, v2, d, U).coeff_monomial(sp.Poly(basis[15], v2, d, U).monoms()[0]))
# No pure kinetic slots either.
pure_kinetic_monomials = [v2**4, v2**3*d**2, v2**2*d**4, v2*d**6, d**8]
for j, mon in enumerate(pure_kinetic_monomials, start=1):
    coeff = sp.Poly(L_dem, v2, d, U).coeff_monomial(sp.Poly(mon, v2, d, U).monoms()[0])
    expect_zero(f"pure kinetic slot l{j}", coeff)

# Axisymmetric variables.
ubar2, a2, b2 = sp.symbols("ubar2 a2 b2", real=True)
subs_axis = {
    c20: ubar2 + 4 * a2,
    c21: ubar2 - a2 + b2,
    c22: ubar2 - a2 - b2,
}
coeffs_axis = {i: sp.simplify(coeffs[i].subs(subs_axis)) for i in coeffs}
print("\nAxisymmetric map:")
for i in range(6, 15):
    print(f"l{i} =", coeffs_axis[i])

print("\nIsotropic specialization a2=b2=0:")
for i in range(6, 15):
    print(f"l{i} =", sp.simplify(coeffs_axis[i].subs({a2: 0, b2: 0})))


banner("PART III — RANK-3 INTERACTION IMAGE AND SIX EXACT SLOT RELATIONS")

M = sp.Matrix(
    [
        [sp.expand(coeffs[i]).coeff(c20), sp.expand(coeffs[i]).coeff(c21), sp.expand(coeffs[i]).coeff(c22)]
        for i in range(6, 15)
    ]
)
print("Interaction map matrix M =")
sp.pprint(M)
print("rank(M) =", M.rank())
if M.rank() != 3:
    raise AssertionError("Unexpected grouped-P2 interaction-map rank.")

left_null = M.T.nullspace()
print("left-nullity =", len(left_null))
if len(left_null) != 6:
    raise AssertionError("Unexpected left-nullity for the 9x3 interaction map.")

L6, L7, L8, L9, L10, L11, L12, L13, L14 = sp.symbols("L6:15")
relations = []
for vec in left_null:
    rel = sp.simplify(sum(vec[i] * [L6, L7, L8, L9, L10, L11, L12, L13, L14][i] for i in range(9)))
    relations.append(rel)

# Canonical readable relation basis.
canonical_relations = [
    2 * L6 + 2 * L7 + L8,
    -L6 - L7 + L9,
    L10 + 2 * L6,
    L11 + L12 - 2 * L6,
    L13 - L6,
    L14 + L11 / 6,
]
print("\nA convenient exact relation basis is:")
for rel in canonical_relations:
    print("  ", rel)

# Verify each canonical relation annihilates the image.
coeff_list = [coeffs[i] for i in range(6, 15)]
subs_image = dict(zip([L6, L7, L8, L9, L10, L11, L12, L13, L14], coeff_list))
for k, rel in enumerate(canonical_relations, start=1):
    expect_zero(f"relation {k} on image", rel.subs(subs_image))


banner("PART IV — EXACT GR 3PN TARGET VIOLATES ALL SIX RELATIONS")

nu = sp.symbols("nu", real=True)
Delta = {
    6: nu * (38 - 116 * nu - 57 * nu**2) / 16,
    7: nu**2 * (20 - 69 * nu) / 16,
    8: 3 * nu**2 * (3 - 11 * nu) / 16,
    9: 5 * nu**3 / 16,
    10: nu * (129 - 98 * nu + 52 * nu**2) / 16,
    11: nu * (-3 + 52 * nu + 124 * nu**2) / 16,
    12: nu * (-5 + 11 * nu + 48 * nu**2) / 12,
    13: -nu * (244 + 3 * sp.pi**2 + 1272 * nu + 96 * nu**2) / 192,
    14: nu * (452 + 3 * sp.pi**2 - 384 * nu - 224 * nu**2) / 64,
}
subs_target = {L6: Delta[6], L7: Delta[7], L8: Delta[8], L9: Delta[9], L10: Delta[10],
               L11: Delta[11], L12: Delta[12], L13: Delta[13], L14: Delta[14]}

violations = []
for k, rel in enumerate(canonical_relations, start=1):
    v = sp.simplify(sp.expand(rel.subs(subs_target)))
    violations.append(v)
    print(f"target violation {k} =", v)
    if v == 0:
        raise AssertionError("A supposed obstruction relation vanished on the exact target.")

nu_eq = sp.Rational(1, 4)
print("\nEqual-mass violations at nu=1/4:")
for k, v in enumerate(violations, start=1):
    print(f"  violation {k} =", sp.simplify(v.subs(nu, nu_eq)))


banner("PART V — FINAL LEDGER")
print("1. The direct grouped-P2 O(w^2) front-end has uniform virial weight 5, so by itself it")
print("   is not a formal 3PN ordinary-Lagrangian block.")
print("2. The unique local isotropic one-step demotion by one inverse orbital weight is 1/U ~ r.")
print("3. After that demotion, the grouped front-end lands exactly in the 9 interaction slots")
print("   l6..l14 with a rank-3 linear image and zero support for l1..l5 or l15.")
print("4. The demoted map obeys six exact slot relations:")
print("      2 l6 + 2 l7 + l8 = 0")
print("      l9 = l6 + l7")
print("      l10 = -2 l6")
print("      l11 + l12 = 2 l6")
print("      l13 = l6")
print("      l14 = -l11/6")
print("5. The exact GR 3PN residual violates all six relations, so the minimal local demoted")
print("   grouped-P2 scaffold cannot be the full 3PN dictionary.")
print("6. Therefore the actual 3PN grouped-P2 constitutive lift must be richer than a single")
print("   local demotion of the direct channel norms; it must introduce additional middle-block")
print("   mixing, and separate mechanisms are still needed for the l1 and l15 sectors.")
