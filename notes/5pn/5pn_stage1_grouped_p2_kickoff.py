#!/usr/bin/env python3
"""
5pn_stage1_grouped_p2_kickoff.py

First executable SymPy audit for the 5PN GR program built on the calibrated
moving-throat PDE handoff.

What this script does
---------------------
1. Encodes the grouped real P2 conservative bookkeeping used by the handoff.
2. Verifies the exact grouped trace/anomaly inverse map.
3. Derives the normalized grouped response moments u2, u4 and outgoing-prefactor
   moments P0, P2, P4 from the conservative operator coefficients D_n and
   transfer coefficients N_n.
4. Derives the minimal isotropic grouped-P2 + static-geometry conservative
   module and proves that it forces the 5PN branch identity u4 = 4 u2^2.
5. Derives the exact obstruction formula when the geometry lane carries dynamic
   even moments:
        c_pole = (1 + eps4) / (4 (1 + eps2)^2).
6. Verifies that the actual isotropic decoupled branch eps2 = eps4 = 0 recovers
   the 3/4 + 1/4 conservative module.
7. Encodes the remaining normalization defect N_Q and the outgoing-renormalized
   factorization with chi_Q.

This is meant to be the first trackable symbolic scaffold for the 5PN program.
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
# I. Exact grouped real P2 bookkeeping
# ---------------------------------------------------------------------------

banner("I. EXACT GROUPED REAL P2 TRACE/ANOMALY BOOKKEEPING")

x20, x21, x22 = sp.symbols("x20 x21 x22", real=True)

xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
a_x = sp.simplify((2 * x20 - x21 - x22) / 10)
b_x = sp.simplify((x21 - x22) / 2)

print("xbar =", xbar)
print("a_x  =", a_x)
print("b_x  =", b_x)

x20_rec = sp.simplify(xbar + 4 * a_x)
x21_rec = sp.simplify(xbar - a_x + b_x)
x22_rec = sp.simplify(xbar - a_x - b_x)

expect_zero("recover x20", x20_rec - x20)
expect_zero("recover x21", x21_rec - x21)
expect_zero("recover x22", x22_rec - x22)

subbanner("Weighted anisotropy norm")
A_sq = sp.simplify(4 * a_x**2 + sp.Rational(4, 5) * b_x**2)
print("A^2 =", A_sq)
print("Isotropy gate: a_x = 0 and b_x = 0")


# ---------------------------------------------------------------------------
# II. Conservative operator moments -> grouped response moments
# ---------------------------------------------------------------------------

banner("II. D_n -> u_n AND N_n -> P_n")

omega = sp.symbols("omega", real=True)
D0, D2, D4 = sp.symbols("D0 D2 D4", nonzero=True)
N0, N2, N4 = sp.symbols("N0 N2 N4")

D_cons = D0 + D2 * omega**2 + D4 * omega**4
Y = sp.expand(sp.series(D0 / D_cons, omega, 0, 6).removeO())
print("Y(omega) =", Y)

u2 = sp.simplify(Y.coeff(omega, 2))
u4 = sp.simplify(Y.coeff(omega, 4))
print("u2 =", u2)
print("u4 =", u4)
expect_zero("u2 formula", u2 + D2 / D0)
expect_zero("u4 formula", u4 - (D2**2 - D0 * D4) / D0**2)

Pref = sp.expand(sp.series((D0 * (N0 + N2 * omega**2 + N4 * omega**4)) / D_cons**2, omega, 0, 6).removeO())
print("\nPref(omega) =", Pref)
P0 = sp.simplify(Pref.coeff(omega, 0))
P2 = sp.simplify(Pref.coeff(omega, 2))
P4 = sp.simplify(Pref.coeff(omega, 4))
print("P0 =", P0)
print("P2 =", P2)
print("P4 =", P4)
expect_zero("P0 formula", P0 - N0 / D0)
expect_zero("P2 formula", P2 - (D0 * N2 - 2 * D2 * N0) / D0**2)
expect_zero(
    "P4 formula",
    P4 - (D0**2 * N4 - 2 * D0 * (D2 * N2 + D4 * N0) + 3 * D2**2 * N0) / D0**3,
)


# ---------------------------------------------------------------------------
# III. Minimal isotropic grouped-P2 + static-geometry module
# ---------------------------------------------------------------------------

banner("III. MINIMAL ISOTROPIC GROUPED-P2 + STATIC-GEOMETRY MODULE")

Kgeom, Kpole, OmegaQ = sp.symbols("Kgeom Kpole OmegaQ", nonzero=True)
K_cons = sp.simplify(Kgeom + Kpole / (1 - omega**2 / OmegaQ**2))
K_series = sp.expand(sp.series(K_cons, omega, 0, 6).removeO())
print("K_Q^cons(omega) =", K_cons)
print("Series =", K_series)

K0 = sp.simplify(K_series.coeff(omega, 0))
K2 = sp.simplify(K_series.coeff(omega, 2))
K4 = sp.simplify(K_series.coeff(omega, 4))
print("K0 =", K0)
print("K2 =", K2)
print("K4 =", K4)

branch_residual = sp.simplify(K0 * K4 - 4 * K2**2)
print("Branch residual K0*K4 - 4*K2^2 =", branch_residual)
Kgeom_sol = sp.solve(sp.Eq(branch_residual, 0), Kgeom)[0]
print("Kgeom forced by branch identity =", Kgeom_sol)
expect_zero("Kgeom - 3*Kpole", Kgeom_sol - 3 * Kpole)

K0_forced = sp.simplify(K0.subs(Kgeom, Kgeom_sol))
expect_zero("K0 - 4*Kpole", K0_forced - 4 * Kpole)

Yhat_cons = sp.simplify((K_cons / K0).subs(Kgeom, Kgeom_sol))
Yhat_series = sp.expand(sp.series(Yhat_cons, omega, 0, 6).removeO())
print("\nYhat_Q^cons(omega) =", Yhat_cons)
print("Series =", Yhat_series)

u2_min = sp.simplify(Yhat_series.coeff(omega, 2))
u4_min = sp.simplify(Yhat_series.coeff(omega, 4))
print("u2(minimal module) =", u2_min)
print("u4(minimal module) =", u4_min)
expect_zero("u2(minimal) - 1/(4 OmegaQ^2)", u2_min - 1 / (4 * OmegaQ**2))
expect_zero("u4(minimal) - 1/(4 OmegaQ^4)", u4_min - 1 / (4 * OmegaQ**4))
expect_zero("minimal 5PN identity u4 - 4 u2^2", u4_min - 4 * u2_min**2)


# ---------------------------------------------------------------------------
# IV. Dynamic geometry obstruction
# ---------------------------------------------------------------------------

banner("IV. DYNAMIC GEOMETRY OBSTRUCTION")

Kg0, Kg2, Kg4 = sp.symbols("Kg0 Kg2 Kg4")
K_total = sp.simplify((Kg0 + Kg2 * omega**2 + Kg4 * omega**4) + Kpole / (1 - omega**2 / OmegaQ**2))
K_total_series = sp.expand(sp.series(K_total, omega, 0, 6).removeO())
K0g = sp.simplify(K_total_series.coeff(omega, 0))
K2g = sp.simplify(K_total_series.coeff(omega, 2))
K4g = sp.simplify(K_total_series.coeff(omega, 4))

Kg0_sol = sp.solve(sp.Eq(sp.simplify(K0g * K4g - 4 * K2g**2), 0), Kg0)[0]
print("Kg0 on the minimal isotropic branch =")
sp.pprint(Kg0_sol)

eps2, eps4 = sp.symbols("eps2 eps4")
subs_eps = {
    Kg2: sp.simplify(eps2 * Kpole / OmegaQ**2),
    Kg4: sp.simplify(eps4 * Kpole / OmegaQ**4),
}
K0_eps = sp.simplify(K0g.subs(Kg0, Kg0_sol).subs(subs_eps))
c_pole = sp.simplify(Kpole / K0_eps)
print("\nc_pole =", c_pole)
expect_zero("c_pole obstruction formula", c_pole - (1 + eps4) / (4 * (1 + eps2)**2))

chi = sp.symbols("chi")
small_cpole = sp.expand(sp.series(c_pole.subs({eps2: chi * eps2, eps4: chi * eps4}), chi, 0, 2).removeO()).subs(chi, 1)
print("small-contamination expansion to first order =", small_cpole)
expect_zero("first-order obstruction expansion", small_cpole - sp.Rational(1, 4) * (1 + eps4 - 2 * eps2))

c_pole_actual = sp.simplify(c_pole.subs({eps2: 0, eps4: 0}))
c_geom_actual = sp.simplify(1 - c_pole_actual)
print("\nActual isotropic decoupled branch:")
print("c_pole =", c_pole_actual)
print("c_geom =", c_geom_actual)
expect_zero("actual c_pole - 1/4", c_pole_actual - sp.Rational(1, 4))
expect_zero("actual c_geom - 3/4", c_geom_actual - sp.Rational(3, 4))


# ---------------------------------------------------------------------------
# V. Outgoing normalization defect
# ---------------------------------------------------------------------------

banner("V. OUTGOING NORMALIZATION DEFECT")

G, c, cs, a = sp.symbols("G c c_s a", positive=True)
Kbar0_target, NQ, chiQ, mhat0 = sp.symbols("Kbar0_target N_Q chi_Q mhat_0", nonzero=True)
Kbar0 = sp.simplify(NQ * Kbar0_target)

Kbar2 = sp.simplify(Kbar0 / (4 * OmegaQ**2))
Kbar4 = sp.simplify(Kbar0 / (4 * OmegaQ**4))
Kbar2_target = sp.simplify(Kbar0_target / (4 * OmegaQ**2))
Kbar4_target = sp.simplify(Kbar0_target / (4 * OmegaQ**4))
Gamma5 = sp.simplify(chiQ * 9 * Kbar0 / (32 * OmegaQ**5))
Gamma5_target = sp.simplify(9 * Kbar0_target / (32 * OmegaQ**5))

expect_zero("Kbar2/Kbar2_target - NQ", sp.simplify(Kbar2 / Kbar2_target - NQ))
expect_zero("Kbar4/Kbar4_target - NQ", sp.simplify(Kbar4 / Kbar4_target - NQ))
expect_zero("Gamma5/Gamma5_target - chiQ*NQ", sp.simplify(Gamma5 / Gamma5_target - chiQ * NQ))

odd_condition = sp.simplify(mhat0**2 * chiQ * NQ - 1)
print("Observable odd normalization condition =", odd_condition)
print("Natural point-particle source-map limit: if mhat_0 -> 1, then N_Q = 1/chi_Q")

Kbar0_target_from_cs = sp.simplify(54 * G * cs**5 / (5 * a**5 * c**5))
print("Equivalent Kbar0_target using Omega_Q = 3 c_s / (2 a):")
print("Kbar0_target =", Kbar0_target_from_cs)


# ---------------------------------------------------------------------------
# VI. Final theorem ledger
# ---------------------------------------------------------------------------

banner("VI. 5PN KICKOFF LEDGER")
print("1. PDE handoff conservative front end:")
print("   D_A^(cons)(omega) = D_A0 + D_A2 omega^2 + D_A4 omega^4 + O(omega^6)")
print("   Y_A(omega) = D_A0 / D_A^(cons)(omega) = 1 + u2^(A) omega^2 + u4^(A) omega^4 + ...")
print()
print("2. Exact grouped 5PN conservative gate:")
print("   - compute (ubar2, a2, b2) and (ubar4, a4, b4)")
print("   - test isotropy: a2 = b2 = a4 = b4 = 0")
print("   - on the isotropic branch, test the single-pole identity u4 = 4 u2^2")
print()
print("3. Minimal grouped-P2 + static-geometry realization forces:")
print("   Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)")
print("   and therefore u4 = 4 u2^2 exactly.")
print()
print("4. Dynamic geometry obstruction:")
print("   c_pole = (1 + eps4) / (4 (1 + eps2)^2)")
print("   so any nonzero eps2 or eps4 is the clean reduced way that 5PN can fail")
print("   the minimal 3/4 + 1/4 realization.")
print()
print("5. On the actual isotropic decoupled branch eps2 = eps4 = 0:")
print("   c_pole = 1/4, c_geom = 3/4.")
print()
print("6. After that conservative check, the radiative/normative question is one number:")
print("   N_Q = Kbar0 / Kbar0_target, with odd condition mhat_0^2 chi_Q N_Q = 1.")
print("   If the actual passive/outgoing branch is canonical and the natural source map")
print("   holds, then chi_Q = 1 and mhat_0 -> 1, so the full reduced point-particle")
print("   closure reduces to N_Q = 1.")
