#!/usr/bin/env python3
"""
2PN inner-throat modal DtN from PDE-side scaffolding (SymPy)
-----------------------------------------------------------
Purpose
-------
This prototype does two things.

1) It builds the simplest genuine PDE-side *unit-test branch*: the regular
   interior Helmholtz Dirichlet-to-Neumann (DtN) operator for a 4D ball.
   That branch is diagonal in harmonic channels and its small-k expansion can
   be written explicitly.

2) It then shows what the minimal *axisymmetric pole completion* has to look
   like in order to reproduce the already-solved conservative 2PN cross sector.

Main takeaways
--------------
A pure isotropic 4D-ball branch is not enough:
- it leaves the monopole channel with zero static stiffness,
- and it cannot split the dipole m=0 and |m|=1 channels.

So the full throat must at least include:
- an O(2)-symmetric (axisymmetric) support sector,
- a regularized monopole support/drive channel,
- and one pure-U geometry closure channel.

The minimal low-frequency completion is therefore the one-pole-per-channel
model with channel admittances
    Y_{1\perp}, Y_{10}, Y_0, Y_{20}, Y_{21}, Y_{22}, Y_g,
whose static residues are already fixed by the solved 2PN derivation.
"""
from __future__ import annotations

import sympy as sp


def banner(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def assert_zero(name: str, expr) -> None:
    res = sp.simplify(sp.expand(expr))
    ok = False
    if isinstance(res, sp.MatrixBase):
        ok = all(sp.simplify(entry) == 0 for entry in res)
    else:
        ok = (res == 0)
    if not ok:
        raise AssertionError(f"{name} failed: {res}")
    print(f"PASS: {name}")


# ---------------------------------------------------------------------------
# 1) PDE unit test: regular interior Helmholtz DtN on a 4D ball
# ---------------------------------------------------------------------------

def ball_dtn_4d_series(ell: int, z, a, order: int = 10):
    """
    Regular interior Helmholtz DtN eigenvalue on the 4D ball for harmonic degree ell.

    In 4 spatial dimensions, the regular radial mode is r^{-1} J_{ell+1}(k r).
    Writing z = k a, the boundary DtN eigenvalue at r = a is
        Lambda_ell(z) = -1/a + (z/a) J'_{ell+1}(z) / J_{ell+1}(z).

    The function returns its exact small-z series.
    """
    nu = ell + 1
    expr = -sp.Integer(1) / a + (z / a) * sp.diff(sp.besselj(nu, z), z) / sp.besselj(nu, z)
    return sp.expand(sp.series(expr, z, 0, order).removeO())


# Symbols
z, a = sp.symbols("z a", positive=True)
omega = sp.symbols("omega", real=True)

Lam0 = ball_dtn_4d_series(0, z, a, 10)
Lam1 = ball_dtn_4d_series(1, z, a, 10)
Lam2 = ball_dtn_4d_series(2, z, a, 10)

Lam0_expected = sp.expand(-z**2 / (4 * a) - z**4 / (96 * a) - z**6 / (1536 * a) - z**8 / (23040 * a))
Lam1_expected = sp.expand(1 / a - z**2 / (6 * a) - z**4 / (288 * a) - z**6 / (8640 * a) - 7 * z**8 / (1658880 * a))
Lam2_expected = sp.expand(2 / a - z**2 / (8 * a) - z**4 / (640 * a) - z**6 / (30720 * a) - 13 * z**8 / (17203200 * a))

banner("1) PDE unit test: isotropic 4D-ball DtN branch")
assert_zero("4D-ball monopole DtN series", Lam0 - Lam0_expected)
assert_zero("4D-ball dipole DtN series", Lam1 - Lam1_expected)
assert_zero("4D-ball quadrupole DtN series", Lam2 - Lam2_expected)
assert_zero("Laplace limit: Lambda_0(0)=0", sp.simplify(Lam0.subs(z, 0)))
assert_zero("Laplace limit: Lambda_1(0)=1/a", sp.simplify(Lam1.subs(z, 0) - 1 / a))
assert_zero("Laplace limit: Lambda_2(0)=2/a", sp.simplify(Lam2.subs(z, 0) - 2 / a))

Y1_ball = sp.expand(sp.series(1 / Lam1, z, 0, 5).removeO())
Y2_ball = sp.expand(sp.series(1 / Lam2, z, 0, 5).removeO())
Y1_ball_expected = sp.expand(a + a * z**2 / 6 + a * z**4 / 32)
Y2_ball_expected = sp.expand(a / 2 + a * z**2 / 32 + 3 * a * z**4 / 1280)
assert_zero("4D-ball dipole admittance series", Y1_ball - Y1_ball_expected)
assert_zero("4D-ball quadrupole admittance series", Y2_ball - Y2_ball_expected)

print("\n4D-ball DtN low-k series:")
print(f"  Lambda_0(z) = {Lam0}")
print(f"  Lambda_1(z) = {Lam1}")
print(f"  Lambda_2(z) = {Lam2}")
print("\nImplications of the isotropic branch:")
print("  * the monopole channel has zero static stiffness (Lambda_0(0)=0)")
print("  * the branch is diagonal in harmonic degree and therefore cannot split dipole m=0 vs |m|=1")
print("  * so a pure isotropic cavity is only a unit test, not the full 2PN throat")


# ---------------------------------------------------------------------------
# 2) Minimal axisymmetric pole completion for the solved 2PN throat operator
# ---------------------------------------------------------------------------

banner("2) Minimal axisymmetric one-pole-per-channel completion")

Om1p, Om10, Om0, Om20, Om21, Om22, Omg = sp.symbols(
    "Omega1p Omega10 Omega0 Omega20 Omega21 Omega22 Omeg", positive=True
)

# Static residues fixed by the solved 2PN program.
R1p = sp.Rational(7, 2)
R10 = sp.Integer(4)
R0 = sp.Integer(1)
R20 = sp.Integer(1)
R21 = sp.Integer(1)
R22 = sp.Integer(1)
Delta = sp.Rational(281, 80)
J0 = sp.Integer(4) / sp.sqrt(5)
J20 = sp.Rational(5, 4)

Y1p = R1p / (1 - omega**2 / Om1p**2)
Y10 = R10 / (1 - omega**2 / Om10**2)
Y0 = R0 / (1 - omega**2 / Om0**2)
Y20 = R20 / (1 - omega**2 / Om20**2)
Y21 = R21 / (1 - omega**2 / Om21**2)
Y22 = R22 / (1 - omega**2 / Om22**2)
Yg = sp.Integer(1) / (1 - omega**2 / Omg**2)

Z1p = sp.simplify(1 / Y1p)
Z10 = sp.simplify(1 / Y10)
Z0 = sp.simplify(1 / Y0)
Z20 = sp.simplify(1 / Y20)
Z21 = sp.simplify(1 / Y21)
Z22 = sp.simplify(1 / Y22)
Zg = sp.simplify(1 / Yg)

assert_zero("Y_{1perp}(0)=7/2", sp.simplify(Y1p.subs(omega, 0) - sp.Rational(7, 2)))
assert_zero("Y_{10}(0)=4", sp.simplify(Y10.subs(omega, 0) - 4))
assert_zero("Y_0(0)=Y_20(0)=Y_21(0)=Y_22(0)=1",
            sp.Matrix([
                sp.simplify(Y0.subs(omega, 0) - 1),
                sp.simplify(Y20.subs(omega, 0) - 1),
                sp.simplify(Y21.subs(omega, 0) - 1),
                sp.simplify(Y22.subs(omega, 0) - 1),
            ]))
assert_zero("Y_g(0)=1", sp.simplify(Yg.subs(omega, 0) - 1))
assert_zero("Z_{1perp}(0)=2/7", sp.simplify(Z1p.subs(omega, 0) - sp.Rational(2, 7)))
assert_zero("Z_{10}(0)=1/4", sp.simplify(Z10.subs(omega, 0) - sp.Rational(1, 4)))
assert_zero("Z_0(0)=Z_20(0)=Z_21(0)=Z_22(0)=1",
            sp.Matrix([
                sp.simplify(Z0.subs(omega, 0) - 1),
                sp.simplify(Z20.subs(omega, 0) - 1),
                sp.simplify(Z21.subs(omega, 0) - 1),
                sp.simplify(Z22.subs(omega, 0) - 1),
            ]))
assert_zero("Z_g(0)=1", sp.simplify(Zg.subs(omega, 0) - 1))

# Even response matrix in the canonical {P0, P20, P21c, P21s, P22c, P22s, U} basis.
M_even = sp.Matrix([
    [Y0, 0, 0, 0, 0, 0, J0 * Y0],
    [0, Y20, 0, 0, 0, 0, J20 * Y20],
    [0, 0, Y21, 0, 0, 0, 0],
    [0, 0, 0, Y21, 0, 0, 0],
    [0, 0, 0, 0, Y22, 0, 0],
    [0, 0, 0, 0, 0, Y22, 0],
    [J0 * Y0, J20 * Y20, 0, 0, 0, 0, J0**2 * Y0 + J20**2 * Y20 - Delta * Yg],
])

R_support = sp.Matrix([
    [sp.sqrt(Y0), 0, 0, 0, 0, 0, J0 * sp.sqrt(Y0)],
    [0, sp.sqrt(Y20), 0, 0, 0, 0, J20 * sp.sqrt(Y20)],
    [0, 0, sp.sqrt(Y21), 0, 0, 0, 0],
    [0, 0, 0, sp.sqrt(Y21), 0, 0, 0],
    [0, 0, 0, 0, sp.sqrt(Y22), 0, 0],
    [0, 0, 0, 0, 0, sp.sqrt(Y22), 0],
])
R_geom = sp.Matrix([[0, 0, 0, 0, 0, 0, sp.sqrt(Delta * Yg)]])
assert_zero("M_even(omega) = support PSD block minus one pure-U closure block",
            R_support.T * R_support - R_geom.T * R_geom - M_even)

M_even_static_expected = sp.Matrix([
    [1, 0, 0, 0, 0, 0, J0],
    [0, 1, 0, 0, 0, 0, J20],
    [0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 1, 0],
    [J0, J20, 0, 0, 0, 0, sp.Rational(5, 4)],
])
assert_zero("Static even matrix reproduces the solved 2PN support/closure data",
            M_even.subs(omega, 0) - M_even_static_expected)

# Effective U-source vector and U-U coefficient.
Jeff = sp.Matrix([J0 * Y0, J20 * Y20, 0, 0, 0, 0])
Seff = sp.simplify(J0**2 * Y0 + J20**2 * Y20 - Delta * Yg)
assert_zero("Static U-U coefficient is 5/4", sp.simplify(Seff.subs(omega, 0) - sp.Rational(5, 4)))

print("\nBare channel DtN kernels:")
print(f"  Z1perp(omega) = {Z1p}")
print(f"  Z10(omega)    = {Z10}")
print(f"  Z0(omega)     = {Z0}")
print(f"  Z20(omega)    = {Z20}")
print(f"  Z21(omega)    = {Z21}")
print(f"  Z22(omega)    = {Z22}")
print(f"  Zg(omega)     = {Zg}")

print("\nFrequency-dependent scalar source data:")
print(f"  Jeff(omega) = {Jeff.T}")
print(f"  Seff(omega) = {Seff}")


# ---------------------------------------------------------------------------
# 3) Static limit of the dynamic model reproduces the solved full cross block
# ---------------------------------------------------------------------------

banner("3) Static limit reproduces the solved conservative 2PN cross operator")

vAx, vAy, vAz, vBx, vBy, vBz = sp.symbols("vAx vAy vAz vBx vBy vBz", real=True)
UA, UB = sp.symbols("UA UB", real=True)

vA2 = sp.expand(vAx**2 + vAy**2 + vAz**2)
vB2 = sp.expand(vBx**2 + vBy**2 + vBz**2)
vAB = sp.expand(vAx * vBx + vAy * vBy + vAz * vBz)
uAB = sp.expand(vAx * vBx + vAy * vBy)
dA = vAz
dB = vBz

Pi0A = sp.expand(sp.sqrt(5) * vA2 / 2)
Pi0B = sp.expand(sp.sqrt(5) * vB2 / 2)
Pi20A = sp.expand((3 * dA**2 - vA2) / 2)
Pi20B = sp.expand((3 * dB**2 - vB2) / 2)
Pi21Ac = sp.expand(sp.sqrt(2) * dA * vAx)
Pi21As = sp.expand(sp.sqrt(2) * dA * vAy)
Pi21Bc = sp.expand(sp.sqrt(2) * dB * vBx)
Pi21Bs = sp.expand(sp.sqrt(2) * dB * vBy)
Pi22Ac = sp.expand((vAx**2 - vAy**2) / (2 * sp.sqrt(2)))
Pi22As = sp.expand((2 * vAx * vAy) / (2 * sp.sqrt(2)))
Pi22Bc = sp.expand((vBx**2 - vBy**2) / (2 * sp.sqrt(2)))
Pi22Bs = sp.expand((2 * vBx * vBy) / (2 * sp.sqrt(2)))

sigma = sp.Rational(1, 2)
eta_perp = sp.Rational(15, 14)
eta_par = sp.Rational(15, 16)

L1wake = sp.expand(-sp.Rational(7, 2) * vAB - sp.Rational(1, 2) * dA * dB)
Lodd_added = sp.expand(
    sigma * (vA2 + vB2) * L1wake
    - (UA + UB) * (eta_perp * sp.Rational(7, 2) * uAB + eta_par * 4 * dA * dB)
)

Leven_dynamic = sp.expand(
    Y0 * (Pi0A + J0 * UA) * (Pi0B + J0 * UB)
    + Y20 * (Pi20A + J20 * UA) * (Pi20B + J20 * UB)
    + Y21 * (Pi21Ac * Pi21Bc + Pi21As * Pi21Bs)
    + Y22 * (Pi22Ac * Pi22Bc + Pi22As * Pi22Bs)
    - Delta * Yg * UA * UB
)

Lfull_static_from_dtn = sp.expand(
    -Y1p * uAB - Y10 * dA * dB + Lodd_added + Leven_dynamic
).subs(omega, 0)

Lfull_static_target = sp.expand(L1wake + Lodd_added + Leven_dynamic.subs(omega, 0))

assert_zero("Static limit of the dynamic model reproduces the solved full cross operator",
            Lfull_static_from_dtn - Lfull_static_target)


# ---------------------------------------------------------------------------
# 4) Low-frequency expansion: what 2PN fixed, what remains PDE data
# ---------------------------------------------------------------------------

banner("4) Low-frequency expansion and remaining PDE observables")

Y_series = {
    "Y1perp": sp.series(Y1p, omega, 0, 4).removeO(),
    "Y10": sp.series(Y10, omega, 0, 4).removeO(),
    "Y0": sp.series(Y0, omega, 0, 4).removeO(),
    "Y20": sp.series(Y20, omega, 0, 4).removeO(),
    "Y21": sp.series(Y21, omega, 0, 4).removeO(),
    "Y22": sp.series(Y22, omega, 0, 4).removeO(),
    "Yg": sp.series(Yg, omega, 0, 4).removeO(),
}
for name, series in Y_series.items():
    print(f"  {name}(omega) = {series}")

Z_series = {
    "Z1perp": sp.series(Z1p, omega, 0, 4).removeO(),
    "Z10": sp.series(Z10, omega, 0, 4).removeO(),
    "Z0": sp.series(Z0, omega, 0, 4).removeO(),
    "Z20": sp.series(Z20, omega, 0, 4).removeO(),
    "Z21": sp.series(Z21, omega, 0, 4).removeO(),
    "Z22": sp.series(Z22, omega, 0, 4).removeO(),
    "Zg": sp.series(Zg, omega, 0, 4).removeO(),
}
for name, series in Z_series.items():
    print(f"  {name}(omega) = {series}")

Seff_series = sp.series(Seff, omega, 0, 4).removeO()
Jeff_series = sp.Matrix([sp.series(entry, omega, 0, 4).removeO() for entry in Jeff])
print("\nLow-frequency source sector:")
print(f"  Jeff(omega) = {Jeff_series.T}")
print(f"  Seff(omega) = {Seff_series}")

Seff_expected = sp.Rational(5, 4) + omega**2 * (
    sp.Rational(16, 5) / Om0**2 + sp.Rational(25, 16) / Om20**2 - sp.Rational(281, 80) / Omg**2
)
assert_zero("Seff(omega) low-frequency coefficient is fixed channel-by-channel",
            Seff_series - Seff_expected)

print("\nDynamic PDE observables still not fixed by conservative 2PN algebra:")
print("  Omega1p, Omega10, Omega0, Omega20, Omega21, Omega22, Omeg")
print("Optional near-spherical reduction:")
print("  set Omega20 = Omega21 = Omega22 to collapse the P2 support poles.")


# ---------------------------------------------------------------------------
# 5) Summary
# ---------------------------------------------------------------------------

banner("5) Summary")
print("1. The isotropic 4D-ball DtN branch is a valid PDE unit test but fails two key checks:")
print("   - no finite static monopole support, and")
print("   - no dipole splitting.")
print("2. The minimal throat completion is therefore axisymmetric and must include:")
print("   - odd dipole channels {1perp, 10},")
print("   - even support channels {0, 20, 21, 22},")
print("   - one pure-U geometry closure channel.")
print("3. The solved 2PN derivation fixes all zero-frequency residues and the scalar source vector,")
print("   but not the seven pole scales. Those are now clean PDE/DtN observables.")
print("\nDone.")
