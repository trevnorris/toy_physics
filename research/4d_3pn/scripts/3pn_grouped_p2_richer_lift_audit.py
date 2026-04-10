#!/usr/bin/env python3
"""
3pn_grouped_p2_richer_lift_audit.py

Audit the next grouped-P2 3PN step after the failed minimal local demotion test.

Main result
-----------
The grouped-P2 ontology itself is not too small.  The failure came from the
minimal constitutive choice.  If one enlarges the local grouped-P2 lift from

    T_A ~ U^{-1} \dot C_A^2

to the natural 9-basis family

    (T_A, S_A, V_A)
    with S_A = U^2 C_A^2,  V_A = (v^2/U) S_A,

then the 9x9 middle-block map into (l6,...,l14) has constant determinant -4/27
and is therefore exactly invertible.  It fits the full solved GR middle block
exactly, while still forcing l1..l5=0 and predicting a specific static
companion l15 = l10+l11+l12+2(l6+l7+l8+l9), which does not equal the true GR
static residual.
"""
from __future__ import annotations
import sympy as sp


def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")


v2, d, U = sp.symbols("v2 d U", real=True)
u2 = v2 - d**2

# Grouped P2 source-square families.
C20sq = sp.expand(sp.Rational(1, 6) * (3 * d**2 - v2 - 2 * U) ** 2)
C21sq = sp.expand(2 * d**2 * u2)
C22sq = sp.expand(sp.Rational(1, 2) * u2**2)

T20 = sp.expand(U * d**2 * (3 * u2 - U) ** 2 / 3)
T21 = sp.expand(U * u2 * (u2 - d**2 - U) ** 2)
T22 = sp.expand(U * d**2 * u2**2)

S20 = sp.expand(U**2 * C20sq)
S21 = sp.expand(U**2 * C21sq)
S22 = sp.expand(U**2 * C22sq)

V20 = sp.expand(v2 * S20 / U)
V21 = sp.expand(v2 * S21 / U)
V22 = sp.expand(v2 * S22 / U)

D20 = sp.expand(d**2 * S20 / U)
D21 = sp.expand(d**2 * S21 / U)
D22 = sp.expand(d**2 * S22 / U)

VT20 = sp.expand(v2 * T20 / U)
VT21 = sp.expand(v2 * T21 / U)
VT22 = sp.expand(v2 * T22 / U)

monoms = {
    1: v2**4,
    2: v2**3 * d**2,
    3: v2**2 * d**4,
    4: v2 * d**6,
    5: d**8,
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


def coeff_vector(expr: sp.Expr) -> list[sp.Expr]:
    poly = sp.Poly(sp.expand(expr), v2, d, U)
    out = []
    for i in range(1, 16):
        mon = sp.Poly(monoms[i], v2, d, U).monoms()[0]
        out.append(sp.simplify(poly.coeff_monomial(mon)))
    return out


families = {
    "T20": T20, "T21": T21, "T22": T22,
    "S20": S20, "S21": S21, "S22": S22,
    "V20": V20, "V21": V21, "V22": V22,
    "D20": D20, "D21": D21, "D22": D22,
    "VT20": VT20, "VT21": VT21, "VT22": VT22,
}


def matrix_for(names: list[str], rows=range(6, 15)) -> sp.Matrix:
    return sp.Matrix([[coeff_vector(families[n])[i - 1] for n in names] for i in rows])


banner("PART I — NATURAL FAMILY RANK TEST")

names_T = ["T20", "T21", "T22"]
names_TS = names_T + ["S20", "S21", "S22"]
names_TSD = names_TS + ["D20", "D21", "D22"]
names_TSV = names_TS + ["V20", "V21", "V22"]
names_TSVT = names_TS + ["VT20", "VT21", "VT22"]

for label, names in [
    ("T", names_T),
    ("T+S", names_TS),
    ("T+S+D", names_TSD),
    ("T+S+VT", names_TSVT),
    ("T+S+V", names_TSV),
]:
    M = matrix_for(names)
    print(f"rank({label}) on (l6..l14) =", M.rank())

A_mid = matrix_for(names_TSV)
det_mid = sp.factor(A_mid.det())
print("det(T+S+V middle-block matrix) =", det_mid)
if det_mid != -sp.Rational(4, 27):
    raise AssertionError("Unexpected determinant for the richer grouped-P2 compiler.")


banner("PART II — EXACT INVERSE MIDDLE-BLOCK COMPILER")

L6, L7, L8, L9, L10, L11, L12, L13, L14 = sp.symbols("L6:15")
generic_target = sp.Matrix([L6, L7, L8, L9, L10, L11, L12, L13, L14])
inv_coeffs = [sp.expand(x) for x in A_mid.LUsolve(generic_target)]

labels = [
    "lambda20", "lambda21", "lambda22",
    "sigma20", "sigma21", "sigma22",
    "tau20", "tau21", "tau22",
]
for label, expr in zip(labels, inv_coeffs):
    print(f"{label} =", expr)

expr_generic = sp.expand(sum(inv_coeffs[i] * families[names_TSV[i]] for i in range(9)))
coords_generic = coeff_vector(expr_generic)

# exact middle-block recovery
for idx, Li in zip(range(6, 15), generic_target):
    expect_zero(f"recovered l{idx} - L{idx}", coords_generic[idx - 1] - Li)

# pure kinetic slots vanish
for idx in range(1, 6):
    expect_zero(f"pure kinetic slot l{idx}", coords_generic[idx - 1])

# static companion
l15_pred = sp.simplify(coords_generic[14])
print("predicted l15 =", l15_pred)
expect_zero(
    "l15 prediction relation",
    l15_pred - (L10 + L11 + L12 + 2 * (L6 + L7 + L8 + L9)),
)


banner("PART III — EXACT FIT TO THE SOLVED GR 3PN MIDDLE BLOCK")

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
    15: nu * (-908 + 63 * sp.pi**2) / 96,
}
target_vec = sp.Matrix([Delta[i] for i in range(6, 15)])
sol_target = [sp.simplify(x) for x in A_mid.LUsolve(target_vec)]
for label, expr in zip(labels, sol_target):
    print(f"{label}(GR) =", expr)

expr_target = sp.expand(sum(sol_target[i] * families[names_TSV[i]] for i in range(9)))
coords_target = coeff_vector(expr_target)
for idx in range(6, 15):
    expect_zero(f"GR middle slot l{idx}", coords_target[idx - 1] - Delta[idx])
for idx in range(1, 6):
    expect_zero(f"GR kinetic slot l{idx}", coords_target[idx - 1])

l15_target_pred = sp.simplify(coords_target[14])
l15_gap = sp.simplify(Delta[15] - l15_target_pred)
print("predicted GR l15 from richer grouped-P2 middle compiler =", l15_target_pred)
print("true GR l15 =", Delta[15])
print("remaining static gap =", l15_gap)


banner("PART IV — TARGET-MINIMALITY WITHIN THE NATURAL (T,S,V) FAMILY")

fit8 = []
for omit in range(len(names_TSV)):
    sub = [names_TSV[i] for i in range(len(names_TSV)) if i != omit]
    A_sub = matrix_for(sub)
    if A_sub.rank() == A_sub.row_join(target_vec).rank():
        fit8.append(sub)

print("8-element fitting subsets =", fit8)
if fit8:
    raise AssertionError("Found an unexpected 8-element subset that still fits the full GR middle block.")


banner("PART V — FINAL LEDGER")
print("1. The demoted dynamic grouped family T_A = U^{-1} dot(C_A)^2 has rank 3 on the")
print("   9 middle slots l6..l14; adjoining the static-support squares S_A = U^2 C_A^2")
print("   lifts this only to rank 6.")
print("2. Among the obvious local scalar dressings by the dimensionless invariants v^2/U")
print("   and d^2/U, the first natural completion that closes the full 9-slot middle block")
print("   is V_A = (v^2/U) S_A = U v^2 C_A^2.")
print("3. The exact 9x9 middle-block compiler built from (T_A,S_A,V_A) has determinant -4/27,")
print("   so it is exactly invertible.")
print("4. Therefore the grouped real P2 ontology can carry the entire solved 3PN middle block")
print("   once this richer local constitutive lift is admitted.")
print("5. The richer grouped compiler still forces l1..l5 = 0 identically, so the pure kinetic")
print("   residual remains outside the grouped-P2 module.")
print("6. The richer grouped compiler predicts a static companion")
print("      l15 = l10 + l11 + l12 + 2(l6+l7+l8+l9),")
print("   which does not equal the true GR static residual; an additional static completion is")
print("   therefore still required.")
