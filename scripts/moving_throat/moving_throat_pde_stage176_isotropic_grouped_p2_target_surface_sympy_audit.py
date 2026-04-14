#!/usr/bin/env python3
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


def expect_zero(name: str, expr) -> None:
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


def grouped_trace_anomaly(x20, x21, x22):
    xbar = sp.simplify((x20 + 2 * x21 + 2 * x22) / 5)
    ax = sp.simplify((2 * x20 - x21 - x22) / 10)
    bx = sp.simplify((x21 - x22) / 2)
    return xbar, ax, bx


def grouped_inverse(xbar, ax, bx):
    x20 = sp.simplify(xbar + 4 * ax)
    x21 = sp.simplify(xbar - ax + bx)
    x22 = sp.simplify(xbar - ax - bx)
    return x20, x21, x22


banner("STAGE 176 — EXACT ISOTROPIC GROUPED-P2 TARGET SURFACE AND THE SCALAR/GEOMETRY FIREWALL")

# ---------------------------------------------------------------------------
# I. Exact grouped conservative trace/anomaly map and isotropic collapse
# ---------------------------------------------------------------------------
subbanner("I. Exact grouped conservative trace/anomaly map")

x20, x21, x22 = sp.symbols("x20 x21 x22", real=True)
xbar, ax, bx = grouped_trace_anomaly(x20, x21, x22)
x20_back, x21_back, x22_back = grouped_inverse(xbar, ax, bx)

print("xbar =")
sp.pprint(xbar)
print("a_x =")
sp.pprint(ax)
print("b_x =")
sp.pprint(bx)
expect_zero("inverse x20", x20_back - x20)
expect_zero("inverse x21", x21_back - x21)
expect_zero("inverse x22", x22_back - x22)

nu2, nu4 = sp.symbols("nu2 nu4", real=True)
a2_iso, b2_iso = grouped_trace_anomaly(nu2, nu2, nu2)[1:]
a4_iso, b4_iso = grouped_trace_anomaly(nu4, nu4, nu4)[1:]
expect_zero("a2 on isotropic common-lane branch", a2_iso)
expect_zero("b2 on isotropic common-lane branch", b2_iso)
expect_zero("a4 on isotropic common-lane branch", a4_iso)
expect_zero("b4 on isotropic common-lane branch", b4_iso)

# ---------------------------------------------------------------------------
# II. Exact one-pole identity in D-space
# ---------------------------------------------------------------------------
subbanner("II. Exact one-pole conservative identity")

D0, D2, D4 = sp.symbols("D0 D2 D4", nonzero=True, real=True)
nu2_common = sp.simplify(-D2 / D0)
nu4_common = sp.simplify((D2**2 - D0 * D4) / D0**2)
Delta_pole = sp.simplify(nu4_common - 4 * nu2_common**2)

print("nu2 =")
sp.pprint(nu2_common)
print("nu4 =")
sp.pprint(nu4_common)
print("Delta_pole =")
sp.pprint(Delta_pole)
expect_zero(
    "Delta_pole + (3 D2^2 + D0 D4)/D0^2",
    Delta_pole + (3 * D2**2 + D0 * D4) / D0**2,
)

D4_onepole = sp.simplify(-3 * D2**2 / D0)
expect_zero(
    "one-pole surface equivalence",
    sp.simplify(Delta_pole.subs(D4, D4_onepole)),
)

# ---------------------------------------------------------------------------
# III. Exact one-parameter carrier on the one-pole surface
# ---------------------------------------------------------------------------
subbanner("III. Exact one-parameter conservative carrier")

omega = sp.symbols("omega", real=True)
OmegaQ2 = sp.simplify(-D0 / (4 * D2))
Y_pole = sp.simplify(sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2 / OmegaQ2))
Y_pole_series = sp.series(Y_pole, omega, 0, 6).removeO()
Y_expected = sp.expand(1 + nu2_common * omega**2 + nu4_common.subs(D4, D4_onepole) * omega**4)

print("Omega_Q^2 =")
sp.pprint(OmegaQ2)
print("Y_pole(omega) =")
sp.pprint(Y_pole)
expect_zero("pole carrier series - expected one-pole series", Y_pole_series - Y_expected)

# ---------------------------------------------------------------------------
# IV. Exact scalar/geometry firewall by Schur complement
# ---------------------------------------------------------------------------
subbanner("IV. Exact scalar/geometry firewall")

K0 = sp.symbols("K0", nonzero=True, real=True)
chi = sp.symbols("chi", real=True)
c20, c21, c22 = sp.symbols("c20 c21 c22", real=True)
Cvec = sp.Matrix([[c20, c21, c22]])
I3 = sp.eye(3)

Deff = sp.simplify(D0 * I3 - chi**2 * Cvec.T * Cvec / K0)
Delta_geom = sp.simplify(Deff - D0 * I3)

print("D_eff,l=2(chi) =")
sp.pprint(Deff)
print("Delta_geom = D_eff - D0 I =")
sp.pprint(Delta_geom)
expect_zero("Delta_geom / chi^2 - expected coefficient", sp.simplify(Delta_geom / chi**2 + (Cvec.T * Cvec) / K0))
expect_zero(
    "d/dchi D_eff at chi=0",
    sp.Matrix(Deff.diff(chi).subs(chi, 0)),
)
expect_zero("D_eff at chi=0 - D0 I", sp.simplify(Deff.subs(chi, 0) - D0 * I3))

banner("STAGE 176 LEDGER")
print("1. The exact grouped trace/anomaly map shows that the isotropic common-lane branch")
print("   is equivalent to a2 = b2 = a4 = b4 = 0 on the conservative grouped bundle.")
print("2. On that isotropic branch the one-pole defect is exactly")
print("      Delta_pole = -(3 D2^2 + D0 D4)/D0^2,"
      " so the one-pole surface is D0 D4 + 3 D2^2 = 0.")
print("3. Therefore the isotropic one-pole conservative carrier is exactly the one-parameter")
print("      3/4 + 1/4 * (1 - omega^2/Omega_Q^2)^(-1)"
      " module through O(omega^4), with Omega_Q^2 = -D0/(4 D2).")
print("4. The scalar/geometry firewall is exact: if l=0 <-> l=2 mixing enters as chi C,")
print("   then the Schur-complement correction to the grouped l=2 block is quadratic,")
print("      D_eff,l=2 = D0 I - chi^2 C^T C / K0,"
      " so there is no O(chi) contamination from the l=0 geometry lane.")
print("5. Stage 176 therefore freezes the first new audited theorem target after Stage 175:")
print("      a2 = b2 = a4 = b4 = 0,   Delta_pole = 0,"
      " together with the exact scalar/geometry firewall.")
