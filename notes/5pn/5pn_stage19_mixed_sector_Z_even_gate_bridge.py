#!/usr/bin/env python3
"""
5pn_stage19_mixed_sector_Z_even_gate_bridge.py

Nineteenth executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Reinstates the conservative Maxwell/mixed Z-sector through the exact invariant
   moments
      Z0 = Q/Delta,
      Z2 = (Q S2 - H Delta)/Delta^2,
      Z4 = [Q(S2^2-Delta) - S2 H Delta]/Delta^3.
2. Derives the exact first-variation formulas
      dZ0, dZ2, dZ4
   in terms of the invariant drifts dQ, dS2, dH, dDelta.
3. Inserts the Stage-13 normalized monomial-rigid kernel and forms the exact
   Z-sector contributions to the conservative even gates
      K1_Z = -dZ2 - dZ0/9,
      H_even,Z = -dZ4 + (2/3) dZ2 + dZ0/27.
4. Evaluates those contributions on a clean positive constructive slice and on
   the pure alpha_W and pure alpha_R directions.
5. Proves that both previously untouched mixed directions are genuinely seen by
   the reinstated Z-sector: the corresponding even-gate values are nonzero.

Interpretation
--------------
Stage 17 left alpha_W and alpha_R untouched only because the Z2,Z4 sector had
been omitted. This script is the first exact reinstatement of that sector. It
shows, on a fully positive constructive branch, that the mixed-sector even gates
become nontrivial precisely along those two directions.
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


def expect_zero(name: str, expr: sp.Expr | sp.Matrix) -> None:
    if isinstance(expr, sp.MatrixBase):
        expr = expr.applyfunc(lambda z: sp.simplify(sp.expand(z)))
        print(f"{name} =")
        sp.pprint(expr)
        if any(entry != 0 for entry in expr):
            raise AssertionError(f"{name} is not zero")
    else:
        expr = sp.simplify(sp.expand(expr))
        print(f"{name} = {expr}")
        if expr != 0:
            raise AssertionError(f"{name} is not zero")


banner("I. EXACT Z-SECTOR INVARIANT BRIDGE")

Q, S2, H, Delta = sp.symbols("Q S2 H Delta", positive=True, real=True)
dQ, dS2, dH, dDelta = sp.symbols("dQ dS2 dH dDelta", real=True)

t = sp.symbols("t", real=True)
Z0_t = (Q + t * dQ) / (Delta + t * dDelta)
Z2_t = ((Q + t * dQ) * (S2 + t * dS2) - (H + t * dH) * (Delta + t * dDelta)) / (Delta + t * dDelta) ** 2
Z4_t = (
    (Q + t * dQ) * ((S2 + t * dS2) ** 2 - (Delta + t * dDelta))
    - (S2 + t * dS2) * (H + t * dH) * (Delta + t * dDelta)
) / (Delta + t * dDelta) ** 3

dZ0 = sp.simplify(sp.diff(Z0_t, t).subs(t, 0))
dZ2 = sp.simplify(sp.diff(Z2_t, t).subs(t, 0))
dZ4 = sp.simplify(sp.diff(Z4_t, t).subs(t, 0))

print("dZ0 =")
sp.pprint(dZ0)
print("\ndZ2 =")
sp.pprint(dZ2)
print("\ndZ4 =")
sp.pprint(dZ4)

expect_zero("dZ0 bridge", dZ0 - (Delta * dQ - Q * dDelta) / Delta**2)
expect_zero(
    "dZ2 bridge",
    dZ2 - (Delta * (-Delta * dH - H * dDelta + Q * dS2 + S2 * dQ) + 2 * dDelta * (Delta * H - Q * S2)) / Delta**3,
)
expect_zero(
    "dZ4 bridge",
    dZ4 + (
        Delta**2 * H * dS2
        + Delta**2 * S2 * dH
        + Delta**2 * dQ
        - 2 * Delta * H * S2 * dDelta
        - 2 * Delta * Q * S2 * dS2
        - 2 * Delta * Q * dDelta
        - Delta * S2**2 * dQ
        + 3 * Q * S2**2 * dDelta
    ) / Delta**4,
)

K1_Z_abstract = sp.simplify(-dZ2 - dZ0 / 9)
He_Z_abstract = sp.simplify(-dZ4 + sp.Rational(2, 3) * dZ2 + dZ0 / 27)

print("\nK1_Z (abstract) =")
sp.pprint(K1_Z_abstract)
print("\nH_even,Z (abstract) =")
sp.pprint(He_Z_abstract)


banner("II. INSERT THE STAGE-13 NORMALIZED KERNEL")

GU, GW, R = sp.symbols("G_U G_W R", positive=True, real=True)
OU, OW = sp.symbols("Omega_U Omega_W", positive=True, real=True)
chi0, deltaU = sp.symbols("chi_0 delta_U", positive=True, real=True)
Estar, Fstar = sp.symbols("E_star F_star", real=True)

aK, aW, aU, aR, aOU = sp.symbols("alpha_K alpha_W alpha_U alpha_R alpha_OmegaU", real=True)

dw = aW
du = aU
dr = aR
dOU = aOU
dDeltaU = sp.simplify(-((1 + deltaU) / (1 + chi0)) * (dr + du - dw - 2 * dOU))
dOW = sp.simplify((dw - du + (1 - Estar) * dOU + Estar * dr - sp.Rational(1, 2) * Fstar * dDeltaU) / (Estar + 2))

Delta_expr = OU**2 * OW**2 - R**2
S2_expr = OU**2 + OW**2
Q_expr = GU**2 * OW**2 + 2 * GU * GW * R + GW**2 * OU**2
H_expr = GU**2 + GW**2

dDelta_expr = sp.simplify(2 * OU**2 * OW**2 * (dOU + dOW) - 2 * R**2 * dr)
dS2_expr = sp.simplify(2 * OU**2 * dOU + 2 * OW**2 * dOW)
dQ_expr = sp.simplify(
    2 * GU**2 * OW**2 * (du + dOW)
    + 2 * GU * GW * R * (du + dw + dr)
    + 2 * GW**2 * OU**2 * (dw + dOU)
)
dH_expr = sp.simplify(2 * GU**2 * du + 2 * GW**2 * dw)

print("Delta =", Delta_expr)
print("S2    =", S2_expr)
print("Q     =", Q_expr)
print("H     =", H_expr)

print("\ndDelta =")
sp.pprint(dDelta_expr)
print("\ndS2 =")
sp.pprint(dS2_expr)
print("\ndQ =")
sp.pprint(dQ_expr)
print("\ndH =")
sp.pprint(dH_expr)

K1_Z = sp.together(K1_Z_abstract.subs({Q: Q_expr, S2: S2_expr, H: H_expr, Delta: Delta_expr, dQ: dQ_expr, dS2: dS2_expr, dH: dH_expr, dDelta: dDelta_expr}))
He_Z = sp.together(He_Z_abstract.subs({Q: Q_expr, S2: S2_expr, H: H_expr, Delta: Delta_expr, dQ: dQ_expr, dS2: dS2_expr, dH: dH_expr, dDelta: dDelta_expr}))

subbanner("Structural point")
print("K1_Z and H_even,Z are exact rational functions of the normalized similarity")
print("variables. They vanish identically only if the conservative Z-sector is omitted.")


banner("III. CONSTRUCTIVE POSITIVE SLICE")

constructive_slice = {
    GU: sp.Integer(5),
    GW: sp.Integer(7),
    R: sp.Integer(2),
    OU: sp.Integer(3),
    OW: sp.Integer(4),
    chi0: sp.Rational(3, 2),
    deltaU: sp.Rational(2, 3),
    Estar: sp.Rational(1, 4),
    Fstar: sp.Rational(5, 6),
}

print("Chosen positive slice:")
for key, val in constructive_slice.items():
    print(f"  {key} = {val}")

Delta_check = sp.simplify(Delta_expr.subs(constructive_slice))
print("\nDelta on the slice =", Delta_check)
if Delta_check <= 0:
    raise AssertionError("The constructive slice must satisfy Delta > 0.")

K1_Z_slice = sp.simplify(K1_Z.subs(constructive_slice))
He_Z_slice = sp.simplify(He_Z.subs(constructive_slice))

print("\nK1_Z on the constructive slice =")
sp.pprint(K1_Z_slice)
print("\nH_even,Z on the constructive slice =")
sp.pprint(He_Z_slice)


banner("IV. THE REINSTATED Z-SECTOR SEES alpha_W AND alpha_R")

pure_aW = {aK: 0, aW: 1, aU: 0, aR: 0, aOU: 0}
pure_aR = {aK: 0, aW: 0, aU: 0, aR: 1, aOU: 0}

K1_aW = sp.simplify(K1_Z_slice.subs(pure_aW))
He_aW = sp.simplify(He_Z_slice.subs(pure_aW))
K1_aR = sp.simplify(K1_Z_slice.subs(pure_aR))
He_aR = sp.simplify(He_Z_slice.subs(pure_aR))

print("Pure alpha_W direction:")
print("  K1_Z     =", K1_aW)
print("  H_even,Z =", He_aW)
print("\nPure alpha_R direction:")
print("  K1_Z     =", K1_aR)
print("  H_even,Z =", He_aR)

if K1_aW == 0 and He_aW == 0:
    raise AssertionError("alpha_W remained invisible after reinstating the Z-sector on the constructive slice.")
if K1_aR == 0 and He_aR == 0:
    raise AssertionError("alpha_R remained invisible after reinstating the Z-sector on the constructive slice.")

print("\nSo the omitted mixed-sector Z2,Z4 moments do exactly what Stage 17 said they")
print("still had to do: they activate the previously untouched alpha_W and alpha_R")
print("directions in the conservative even-gate system.")


banner("V. FINAL LEDGER")
print("1. The conservative Maxwell/mixed moments Z0,Z2,Z4 have exact first-variation")
print("   formulas in the invariant variables (Q,S2,H,Delta).")
print("2. Their contributions to the 5PN even gates are")
print("      K1_Z     = -dZ2 - dZ0/9,")
print("      H_even,Z = -dZ4 + (2/3) dZ2 + dZ0/27.")
print("3. Inserting the Stage-13 normalized monomial-rigid kernel makes both gates")
print("   honest functions of the mixed-sector similarity directions alpha_W and alpha_R.")
print("4. On a clean positive constructive branch, the pure alpha_W and pure alpha_R")
print("   directions give nonzero even-gate values.")
print("5. Therefore the omission in Stage 17 was exactly where the mixed-sector freedom")
print("   was hiding: once Z2,Z4 are reinstated, alpha_W and alpha_R stop being free.")
