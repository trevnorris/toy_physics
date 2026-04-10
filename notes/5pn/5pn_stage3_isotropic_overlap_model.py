#!/usr/bin/env python3
"""
5pn_stage3_isotropic_overlap_model.py

Third executable SymPy audit for the 5PN grouped-real-P2 program.

What this script does
---------------------
1. Chooses one explicit finite-throat radial/axial profile family on s in [0,L]:
      - wall profile chi_eta(s) = sqrt(2/L) sin(pi s / L),
      - D/N-like half-wave profile phi_DN(s) = sqrt(2/L) sin(pi s / (2L)).
2. Computes the exact overlap integrals that feed the isotropic grouped-real-P2
   bundle on the Stage-7 reduced branch.
3. Builds one explicit isotropic one-BdG + one Maxwell/mixed-pair prototype for
   the low-frequency coefficients
      B0,B2,B4, Z0,Z2,Z4, N0,N2,N4.
4. Constructs the conservative operator moments D0,D2,D4 and the normalized
   grouped response / outgoing prefactor moments u2,u4 and P0,P2,P4.
5. Verifies the exact grouped isotropy collapse when all three grouped lanes
   share the same radial/axial data.
6. Derives two sharp algebraic compatibility surfaces:
      (a) the 5PN conservative one-pole identity u4 = 4 u2^2,
      (b) the universal 2.5PN/4PN normalization target for P0.
7. Shows that in this explicit prototype both conditions fix the same static
   wall stiffness K, so simultaneous success reduces to one compatibility
   equation between the overlap bundles.

Interpretation
--------------
This is the first explicit radial/axial overlap model in the 5PN stack. It does
not claim to be the true moving-throat branch; it is the first concrete place
where the handoff's B/Z/N bundle can be evaluated exactly and compared against
both the 5PN conservative identity and the 2.5PN/4PN normalization target.
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


# ---------------------------------------------------------------------------
# I. Explicit finite-throat profile family
# ---------------------------------------------------------------------------

banner("I. EXPLICIT FINITE-THROAT PROFILES AND OVERLAP INTEGRALS")

s, L = sp.symbols("s L", positive=True, real=True)

chi_eta = sp.sqrt(2 / L) * sp.sin(sp.pi * s / L)
phi_DN = sp.sqrt(2 / L) * sp.sin(sp.pi * s / (2 * L))
# Use the wall profile itself as the simplest brane-like U profile, and the
# D/N half-wave as the simplest mixed W profile.
u_prof = chi_eta
w_prof = phi_DN

print("chi_eta(s) =", chi_eta)
print("phi_DN(s)  =", phi_DN)

norm_eta = sp.simplify(sp.integrate(chi_eta**2, (s, 0, L)))
norm_dn = sp.simplify(sp.integrate(phi_DN**2, (s, 0, L)))
I_eta_phi_exact = sp.simplify(sp.integrate(chi_eta * phi_DN, (s, 0, L)))
I_eta_u_exact = sp.simplify(sp.integrate(chi_eta * u_prof, (s, 0, L)))
I_eta_w_exact = sp.simplify(sp.integrate(chi_eta * w_prof, (s, 0, L)))
I_uw_exact = sp.simplify(sp.integrate(u_prof * w_prof, (s, 0, L)))

expect_zero("||chi_eta||^2 - 1", norm_eta - 1)
expect_zero("||phi_DN||^2 - 1", norm_dn - 1)
expect_zero("I_eta_u - 1", I_eta_u_exact - 1)
expect_zero("I_eta_w - I_eta_phi", I_eta_w_exact - I_eta_phi_exact)
expect_zero("I_uw - I_eta_phi", I_uw_exact - I_eta_phi_exact)
expect_zero("I_eta_phi - 8/(3 pi)", I_eta_phi_exact - sp.Rational(8, 3) / sp.pi)

# Keep a symbolic alias for readability in the rest of the script.
I_mix = sp.symbols("I_mix", positive=True, real=True)
print("\nFor the rest of the audit we freeze the exact overlap alias")
print("  I_mix =", I_eta_phi_exact)


# ---------------------------------------------------------------------------
# II. One-BdG + one Maxwell/mixed-pair isotropic prototype
# ---------------------------------------------------------------------------

banner("II. ISOTROPIC ONE-BDG + ONE MAXWELL/MIXED-PAIR PROTOTYPE")

lamB, lamU, lamW, lamR = sp.symbols("lambda_B lambda_U lambda_W lambda_R", real=True)
varpi = sp.symbols("varpi", positive=True, real=True)
OmegaU, OmegaW = sp.symbols("Omega_U Omega_W", positive=True, real=True)
K, M = sp.symbols("K M", positive=True, real=True)

# Effective couplings from the explicit overlap family.
C = lamB * I_mix
GU = lamU
GW = lamW * I_mix
R = lamR * I_mix

print("C  = lambda_B * I_mix")
print("GU = lambda_U")
print("GW = lambda_W * I_mix")
print("R  = lambda_R * I_mix")

# BdG support moments.
B0 = C**2 / varpi**2
B2 = C**2 / varpi**4
B4 = C**2 / varpi**6

# One-lane conservative Maxwell/mixed invariants.
Delta = OmegaU**2 * OmegaW**2 - R**2
S2 = OmegaU**2 + OmegaW**2
Q = GU**2 * OmegaW**2 + 2 * GU * GW * R + GW**2 * OmegaU**2
H = GU**2 + GW**2
Pproto = OmegaU**2 * GW + R * GU

Z0 = Q / Delta
Z2 = (Q * S2 - H * Delta) / Delta**2
Z4 = (Q * (S2**2 - Delta) - S2 * H * Delta) / Delta**3

N0 = Pproto**2 / Delta**2
N2 = 2 * Pproto * (Pproto * S2 - Delta * GW) / Delta**3
N4 = (
    Delta**2 * GW**2
    - 2 * Delta * Pproto**2
    - 4 * Delta * Pproto * S2 * GW
    + 3 * Pproto**2 * S2**2
) / Delta**4

print("\nBdG moments:")
print("B0 =", B0)
print("B2 =", B2)
print("B4 =", B4)

print("\nMaxwell/mixed invariants:")
print("Delta  =", Delta)
print("S2     =", S2)
print("Q      =", Q)
print("H      =", H)
print("Pproto =", Pproto)

print("\nMaxwell/mixed low-frequency moments:")
print("Z0 =", Z0)
print("Z2 =", Z2)
print("Z4 =", Z4)
print("N0 =", N0)
print("N2 =", N2)
print("N4 =", N4)


# ---------------------------------------------------------------------------
# III. Conservative operator and grouped response / prefactor moments
# ---------------------------------------------------------------------------

banner("III. D0,D2,D4 -> u2,u4 AND N0,N2,N4 -> P0,P2,P4")

D0 = K - B0 - Z0
D2 = -(M + B2 + Z2)
D4 = -(B4 + Z4)

u2 = -D2 / D0
u4 = (D2**2 - D0 * D4) / D0**2

P0 = N0 / D0
P2 = (D0 * N2 - 2 * D2 * N0) / D0**2
P4 = (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3

print("D0 =", D0)
print("D2 =", D2)
print("D4 =", D4)

print("\nu2 =", u2)
print("u4 =", u4)
print("P0 =", P0)
print("P2 =", P2)
print("P4 =", P4)

subbanner("Exact grouped isotropy collapse")

u2_20, u2_21, u2_22 = u2, u2, u2
u4_20, u4_21, u4_22 = u4, u4, u4
P0_20, P0_21, P0_22 = P0, P0, P0

ubar2 = (u2_20 + 2 * u2_21 + 2 * u2_22) / 5
a2 = (2 * u2_20 - u2_21 - u2_22) / 10
b2 = (u2_21 - u2_22) / 2

ubar4 = (u4_20 + 2 * u4_21 + 2 * u4_22) / 5
a4 = (2 * u4_20 - u4_21 - u4_22) / 10
b4 = (u4_21 - u4_22) / 2

P0bar = (P0_20 + 2 * P0_21 + 2 * P0_22) / 5
aP0 = (2 * P0_20 - P0_21 - P0_22) / 10
bP0 = (P0_21 - P0_22) / 2

expect_zero("ubar2 - u2", ubar2 - u2)
expect_zero("a2", a2)
expect_zero("b2", b2)
expect_zero("ubar4 - u4", ubar4 - u4)
expect_zero("a4", a4)
expect_zero("b4", b4)
expect_zero("P0bar - P0", P0bar - P0)
expect_zero("aP0", aP0)
expect_zero("bP0", bP0)


# ---------------------------------------------------------------------------
# IV. Two sharp compatibility surfaces
# ---------------------------------------------------------------------------

banner("IV. 5PN ONE-POLE SURFACE AND 2.5PN/4PN NORMALIZATION SURFACE")

# Use compact abbreviations for the positive quantities entering the one-pole test.
A2 = M + B2 + Z2
A4 = B4 + Z4

residual_5pn = (D0 * A4 - 3 * A2**2) / D0**2
print("5PN conservative residual u4 - 4 u2^2 =")
sp.pprint(residual_5pn)

K_from_5pn = B0 + Z0 + 3 * A2**2 / A4
print("\nStatic wall stiffness required by u4 = 4 u2^2:")
print("K_(5PN) =", K_from_5pn)

G, c, c_s, a_th, mhat0 = sp.symbols("G c c_s a_th mhat_0", positive=True, real=True)
P0_target = 54 * G * c_s**5 / (5 * a_th**5 * c**5 * mhat0**2)
K_from_norm = B0 + Z0 + N0 / P0_target

print("\nUniversal normalization target:")
print("P0_target =", P0_target)
print("Static wall stiffness required by the 2.5PN/4PN target:")
print("K_(norm) =", K_from_norm)

compatibility = N0 / P0_target - 3 * A2**2 / A4
print("\nSimultaneous compatibility indicator K_(norm) - K_(5PN) =")
sp.pprint(compatibility)
print("This must vanish if one static K is to satisfy both the 5PN conservative")
print("identity and the 2.5PN/4PN normalization target in this explicit prototype.")

subbanner("Useful equivalent form of the compatibility condition")
compatibility_target = sp.simplify(sp.Eq(N0 / P0_target, 3 * A2**2 / A4))
print(compatibility_target)


# ---------------------------------------------------------------------------
# V. Final ledger
# ---------------------------------------------------------------------------

banner("V. EXPLICIT ISOTROPIC OVERLAP LEDGER")
print("1. Exact finite-throat overlaps on the chosen profile family:")
print("   I_eta_u = 1,")
print("   I_eta_phi = I_eta_w = I_uw = 8/(3 pi).")
print("2. One explicit isotropic prototype then yields closed formulas for")
print("   (B0,B2,B4), (Z0,Z2,Z4), and (N0,N2,N4).")
print("3. The grouped real {20,21,22} bundle collapses exactly:")
print("      a2 = b2 = a4 = b4 = aP0 = bP0 = 0.")
print("4. In this prototype the 5PN conservative identity u4 = 4 u2^2 fixes one")
print("   static wall stiffness combination K_(5PN).")
print("5. The universal 2.5PN/4PN normalization target fixes another one, K_(norm).")
print("6. Therefore simultaneous success is equivalent to one compatibility equation:")
print("      N0 / P0_target = 3 (M + B2 + Z2)^2 / (B4 + Z4).")
print("7. The next honest step is to perturb this explicit isotropic bundle by a weak")
print("   axisymmetric Y20 deformation and compute the grouped obstruction slopes.")
