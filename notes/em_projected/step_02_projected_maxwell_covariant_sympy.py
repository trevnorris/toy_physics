#!/usr/bin/env python3
"""Derive the exact projected inhomogeneous Maxwell law from the localized 4+1 bulk equation.

This script is intentionally symbolic and notation-heavy. It keeps the projection
kernel W(w) explicit and treats the gauge term generically as Gamma^nu, so the
result applies both to the current unweighted gauge-fixing form and to later
weighted variants.

Outputs:
- projection identities
- projected inhomogeneous law for a generic brane component nu
- boundary-decay simplification
- projected charge continuity (open-system charge law)
"""
from __future__ import annotations

import sympy as sp


def line(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.simplify(expr)
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    value = sp.simplify(expr)
    if value == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


# Coordinates and constants

t, x, y, z, w = sp.symbols("t x y z w", real=True)
mu0, xi = sp.symbols("mu0 xi", nonzero=True)
coords_4 = (t, x, y, z)

# Projection kernel and generic bulk fields
W = sp.Function("W")(w)
Wp = sp.diff(W, w)

Q = sp.Function("Q")(t, x, y, z, w)
G0 = sp.Function("G0_nu")(t, x, y, z, w)   # Z F^{0 nu}
Gx = sp.Function("Gx_nu")(t, x, y, z, w)   # Z F^{x nu}
Gy = sp.Function("Gy_nu")(t, x, y, z, w)   # Z F^{y nu}
Gz = sp.Function("Gz_nu")(t, x, y, z, w)   # Z F^{z nu}
Gw = sp.Function("Gw_nu")(t, x, y, z, w)   # Z F^{w nu}
Gamma = sp.Function("Gamma_nu")(t, x, y, z, w)  # (1/xi) partial^nu B  or weighted variant
Jnu = sp.Function("J_nu")(t, x, y, z, w)
J0 = sp.Function("J0")(t, x, y, z, w)
Jx = sp.Function("Jx")(t, x, y, z, w)
Jy = sp.Function("Jy")(t, x, y, z, w)
Jz = sp.Function("Jz")(t, x, y, z, w)
Jw = sp.Function("Jw")(t, x, y, z, w)

ww_int = (w, -sp.oo, sp.oo)


def P(expr: sp.Expr) -> sp.Integral:
    """Project with kernel W(w)."""
    return sp.Integral(W * expr, ww_int)


def Pwprime(expr: sp.Expr) -> sp.Integral:
    """Projection weighted by W'(w), which appears after integrating by parts."""
    return sp.Integral(Wp * expr, ww_int)


def boundary(expr: sp.Expr) -> sp.Integral:
    """Formal total-derivative boundary term in w."""
    return sp.Integral(sp.diff(expr, w), ww_int)


def boundary_value(expr: sp.Expr) -> sp.Expr:
    """Evaluate the actual boundary contribution for a concrete profile."""
    return sp.simplify(sp.limit(expr, w, sp.oo) - sp.limit(expr, w, -sp.oo))


# -----------------------------------------------------------------------------
# Projection identities
# -----------------------------------------------------------------------------
line("1) Projection identities")

print("Projection commutes with t,x,y,z derivatives because W depends only on w:")
print("   Proj_W[partial_t Q] = partial_t Proj_W[Q]")
print("   Proj_W[partial_x Q] = partial_x Proj_W[Q]")
print("   Proj_W[partial_y Q] = partial_y Proj_W[Q]")
print("   Proj_W[partial_z Q] = partial_z Proj_W[Q]")

proj_dw = sp.Eq(P(sp.diff(Q, w)), boundary(W * Q) - Pwprime(Q))
print("\nProjection does NOT commute with the transverse derivative:")
print("  ", sp.sstr(proj_dw))
print("\nIf boundary terms vanish at |w| -> infinity, this reduces to:")
proj_dw_decay = sp.Eq(P(sp.diff(Q, w)), -Pwprime(Q))
print("  ", sp.sstr(proj_dw_decay))


# -----------------------------------------------------------------------------
# Generic projected inhomogeneous equation
# -----------------------------------------------------------------------------
line("2) Exact projected inhomogeneous Maxwell law for a generic brane component nu")

bulk_eq = sp.Eq(
    sp.diff(G0, t) + sp.diff(Gx, x) + sp.diff(Gy, y) + sp.diff(Gz, z) + sp.diff(Gw, w) + Gamma,
    mu0 * Jnu,
)
print("Bulk equation written in generic projected-brane notation:")
print("  ", sp.sstr(bulk_eq))

proj_bulk_eq = sp.Eq(P(bulk_eq.lhs), P(bulk_eq.rhs))
print("\nProject with W(w):")
print("  ", sp.sstr(proj_bulk_eq))

projected_eq = sp.Eq(
    sp.diff(P(G0), t)
    + sp.diff(P(Gx), x)
    + sp.diff(P(Gy), y)
    + sp.diff(P(Gz), z)
    + boundary(W * Gw)
    - Pwprime(Gw)
    + P(Gamma),
    mu0 * P(Jnu),
)
print("\nAfter commuting brane derivatives and integrating the w-derivative by parts:")
print("  ", sp.sstr(projected_eq))

projected_eq_decay = sp.Eq(
    sp.diff(P(G0), t)
    + sp.diff(P(Gx), x)
    + sp.diff(P(Gy), y)
    + sp.diff(P(Gz), z)
    - Pwprime(Gw)
    + P(Gamma),
    mu0 * P(Jnu),
)
print("\nUnder boundary decay / compact support in w:")
print("  ", sp.sstr(projected_eq_decay))

print("\nDictionary back to the localized Maxwell theory:")
print("  G^M_nu  := Z(w) F^{M nu}")
print("  Gamma_nu := (1/xi) partial^nu (partial · A)  [or weighted-gauge variant]")
print("  J_nu     := bulk source current component J^nu")

print("\nSo the exact projected brane equation is:")
print("  partial_mu Proj_W[ Z F^{mu nu} ]")
print("  + Boundary[ W Z F^{w nu} ] - Proj_{W'}[ Z F^{w nu} ]")
print("  + Proj_W[ (1/xi) partial^nu (partial · A) ]")
print("  = mu0 Proj_W[ J^nu ]")

print("\nThis is the precise place where projection differs from reduction:")
print("  - the brane sees an explicit transverse leakage term from Z F^{w nu},")
print("  - the source-coupled flux is Proj_W[ Z F^{mu nu} ], not simply Proj_W[F^{mu nu}],")
print("  - and the gauge-fixing contribution survives as a projected gauge-driver term.")


# -----------------------------------------------------------------------------
# Projected charge continuity
# -----------------------------------------------------------------------------
line("3) Projected charge continuity (open-system charge law)")

bulk_cont = sp.Eq(sp.diff(J0, t) + sp.diff(Jx, x) + sp.diff(Jy, y) + sp.diff(Jz, z) + sp.diff(Jw, w), 0)
print("Bulk current conservation:")
print("  ", sp.sstr(bulk_cont))

proj_cont = sp.Eq(
    sp.diff(P(J0), t)
    + sp.diff(P(Jx), x)
    + sp.diff(P(Jy), y)
    + sp.diff(P(Jz), z)
    + boundary(W * Jw)
    - Pwprime(Jw),
    0,
)
print("\nProjected exact continuity law:")
print("  ", sp.sstr(proj_cont))

proj_cont_decay = sp.Eq(
    sp.diff(P(J0), t)
    + sp.diff(P(Jx), x)
    + sp.diff(P(Jy), y)
    + sp.diff(P(Jz), z),
    Pwprime(Jw),
)
print("\nIf the boundary term vanishes, the projected charge law becomes:")
print("  ", sp.sstr(proj_cont_decay))
print("\nInterpretation: apparent brane charge nonconservation is generated by transverse current J^w.")


# -----------------------------------------------------------------------------
# Compact summary block
# -----------------------------------------------------------------------------
line("4) Compact summary")
print("Exact projected inhomogeneous law (boundary-decay form):")
print("  partial_mu Proj_W[Z F^{mu nu}] - Proj_{W'}[Z F^{w nu}] + Proj_W[Gamma^nu] = mu0 Proj_W[J^nu]")
print("where Gamma^nu is the gauge term and Proj_{W'}[·] := int W'(w) (·) dw.")
print()
print("Exact projected charge continuity (boundary-decay form):")
print("  partial_t Proj_W[J^0] + partial_a Proj_W[J^a] = Proj_{W'}[J^w]")
print()
print("So a projected Maxwell theory is naturally an open-system electrodynamics, not a closed copy of 3+1 Maxwell.")


# -----------------------------------------------------------------------------
# Concrete audit block
# -----------------------------------------------------------------------------
line("5) Concrete symbolic audit")

Wg = sp.exp(-w**2) / sp.sqrt(sp.pi)
Wgp = sp.diff(Wg, w)

def Pg(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.integrate(Wg * expr, (w, -sp.oo, sp.oo)))


def Pgp(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.integrate(Wgp * expr, (w, -sp.oo, sp.oo)))


Phi = (t**2 + x) * (w**2 + 1)
assert_zero("projection commutes with brane derivative", Pg(sp.diff(Phi, t)) - sp.diff(Pg(Phi), t))

Q_concrete = w**3 + w
boundary_Q = boundary_value(Wg * Q_concrete)
assert_zero("decaying-kernel boundary term", boundary_Q)
assert_zero("decaying-kernel integration by parts with boundary", Pg(sp.diff(Q_concrete, w)) - (boundary_Q - Pgp(Q_concrete)))
assert_zero("decaying-kernel integration by parts", Pg(sp.diff(Q_concrete, w)) + Pgp(Q_concrete))
assert_nonzero("mutated decaying-kernel IBP sign should fail", Pg(sp.diff(Q_concrete, w)) - (boundary_Q + Pgp(Q_concrete)))

G0c = t * (w**2 + 1)
Gxc = x * w**2
Gyc = y * (w**2 + 2)
Gzc = z * w**2
Gwc = w
Gammac = (t + x) * (w**2 + 1)

project_bulk_lhs = Pg(
    sp.diff(G0c, t)
    + sp.diff(Gxc, x)
    + sp.diff(Gyc, y)
    + sp.diff(Gzc, z)
    + sp.diff(Gwc, w)
    + Gammac
)
projected_lhs = (
    sp.diff(Pg(G0c), t)
    + sp.diff(Pg(Gxc), x)
    + sp.diff(Pg(Gyc), y)
    + sp.diff(Pg(Gzc), z)
    - Pgp(Gwc)
    + Pg(Gammac)
)
boundary_Gw = boundary_value(Wg * Gwc)
projected_lhs_with_boundary = projected_lhs + boundary_Gw
assert_zero("concrete projected inhomogeneous boundary term", boundary_Gw)
assert_zero("concrete projected inhomogeneous law with boundary", project_bulk_lhs - projected_lhs_with_boundary)
assert_zero("concrete projected inhomogeneous law", project_bulk_lhs - projected_lhs)
assert_nonzero(
    "mutated projected inhomogeneous leakage sign should fail",
    project_bulk_lhs
    - (
        sp.diff(Pg(G0c), t)
        + sp.diff(Pg(Gxc), x)
        + sp.diff(Pg(Gyc), y)
        + sp.diff(Pg(Gzc), z)
        + boundary_Gw
        + Pgp(Gwc)
        + Pg(Gammac)
    ),
)

leakage = -Pgp(Gwc)
assert_zero("concrete transverse leakage value", leakage - 1)
assert_nonzero("concrete transverse leakage", leakage)

J0c = t * (w**2 + 1)
Jxc = x * w**2
Jyc = y * (w**2 + 2)
Jzc = z * w**2
Jwc = w
project_bulk_cont = Pg(
    sp.diff(J0c, t)
    + sp.diff(Jxc, x)
    + sp.diff(Jyc, y)
    + sp.diff(Jzc, z)
    + sp.diff(Jwc, w)
)
projected_cont = (
    sp.diff(Pg(J0c), t)
    + sp.diff(Pg(Jxc), x)
    + sp.diff(Pg(Jyc), y)
    + sp.diff(Pg(Jzc), z)
    - Pgp(Jwc)
)
boundary_Jw = boundary_value(Wg * Jwc)
assert_zero("concrete projected continuity boundary term", boundary_Jw)
assert_zero("concrete projected continuity law with boundary", project_bulk_cont - (projected_cont + boundary_Jw))
assert_zero("concrete projected continuity law", project_bulk_cont - projected_cont)

print("Concrete Gaussian-kernel checks: commutation, explicit boundary discharge, nonzero leakage, and projected continuity all pass.")
print("STATUS: PASS")
