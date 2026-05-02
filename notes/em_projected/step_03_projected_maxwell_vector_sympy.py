#!/usr/bin/env python3
"""Derive projected Maxwell-style vector equations from the projected 4+1 system.

This script separates two kinds of brane objects:
- measured fields from projected F_{MN}
- source-coupled flux fields from projected Z F^{MN}

It then verifies:
1. homogeneous laws project cleanly,
2. inhomogeneous laws acquire transverse leakage and gauge-driver terms,
3. the projected theory naturally distinguishes (E,B) from (D,H).
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


# Brane coordinates

t, x, y, z = sp.symbols("t x y z", real=True)
mu0 = sp.symbols("mu0", nonzero=True)

# Measured fields from projected F_{MN}
E1 = sp.Function("E1_meas")(t, x, y, z)
E2 = sp.Function("E2_meas")(t, x, y, z)
E3 = sp.Function("E3_meas")(t, x, y, z)
B1 = sp.Function("B1_meas")(t, x, y, z)
B2 = sp.Function("B2_meas")(t, x, y, z)
B3 = sp.Function("B3_meas")(t, x, y, z)

# Source-coupled flux fields from projected Z F^{MN}
D1 = sp.Function("D1_flux")(t, x, y, z)
D2 = sp.Function("D2_flux")(t, x, y, z)
D3 = sp.Function("D3_flux")(t, x, y, z)
H1 = sp.Function("H1_flux")(t, x, y, z)
H2 = sp.Function("H2_flux")(t, x, y, z)
H3 = sp.Function("H3_flux")(t, x, y, z)

# Projected source, leakage, and gauge-driver terms
rho = sp.Function("rho_proj")(t, x, y, z)
J1 = sp.Function("J1_proj")(t, x, y, z)
J2 = sp.Function("J2_proj")(t, x, y, z)
J3 = sp.Function("J3_proj")(t, x, y, z)
L0 = sp.Function("Leak0")(t, x, y, z)
L1 = sp.Function("Leak1")(t, x, y, z)
L2 = sp.Function("Leak2")(t, x, y, z)
L3 = sp.Function("Leak3")(t, x, y, z)
G0 = sp.Function("Gauge0")(t, x, y, z)
G1 = sp.Function("Gauge1")(t, x, y, z)
G2 = sp.Function("Gauge2")(t, x, y, z)
G3 = sp.Function("Gauge3")(t, x, y, z)

# -----------------------------------------------------------------------------
# Homogeneous sector from projected Bianchi identities
# -----------------------------------------------------------------------------
line("1) Homogeneous projected equations")

# Definitions consistent with E_a = \bar F_{a0}, B = (\bar F_{23}, \bar F_{31}, \bar F_{12})
F10, F20, F30 = E1, E2, E3
F23, F31, F12 = B1, B2, B3
F01, F02, F03 = -E1, -E2, -E3

# div B = 0 from indices (1,2,3)
divB = sp.diff(F23, x) + sp.diff(F31, y) + sp.diff(F12, z)
print("Measured-field definitions:")
print("  E = ( Proj_W[F_10], Proj_W[F_20], Proj_W[F_30] )")
print("  B = ( Proj_W[F_23], Proj_W[F_31], Proj_W[F_12] )")
print("\nDivergence law from the projected Bianchi identity with indices (1,2,3):")
print("  ", sp.sstr(sp.Eq(divB, 0)))

# Faraday from (2,3,0), (3,1,0), (1,2,0)
far1 = sp.diff(F23, t) + sp.diff(F30, y) + sp.diff(F02, z)
far2 = sp.diff(F31, t) + sp.diff(F10, z) + sp.diff(F03, x)
far3 = sp.diff(F12, t) + sp.diff(F20, x) + sp.diff(F01, y)

curlE1 = sp.diff(E3, y) - sp.diff(E2, z)
curlE2 = sp.diff(E1, z) - sp.diff(E3, x)
curlE3 = sp.diff(E2, x) - sp.diff(E1, y)

checks_faraday = [
    sp.simplify(far1 - (sp.diff(B1, t) + curlE1)),
    sp.simplify(far2 - (sp.diff(B2, t) + curlE2)),
    sp.simplify(far3 - (sp.diff(B3, t) + curlE3)),
]

print("\nFaraday components from the projected Bianchi identities:")
print("  ", sp.sstr(sp.Eq(sp.diff(B1, t) + curlE1, 0)))
print("  ", sp.sstr(sp.Eq(sp.diff(B2, t) + curlE2, 0)))
print("  ", sp.sstr(sp.Eq(sp.diff(B3, t) + curlE3, 0)))
print("\nVerification residues (all should be 0):")
print("  ", checks_faraday)
for i, residue in enumerate(checks_faraday, start=1):
    assert_zero(f"Faraday component {i}", residue)

print("\nCompact homogeneous vector form:")
print("  div B_meas = 0")
print("  curl E_meas + partial_t B_meas = 0")


# -----------------------------------------------------------------------------
# Inhomogeneous sector in vector form
# -----------------------------------------------------------------------------
line("2) Inhomogeneous projected equations")

print("Source-coupled flux-field definitions:")
print("  D = ( Proj_W[ Z F^{10} ], Proj_W[ Z F^{20} ], Proj_W[ Z F^{30} ] )")
print("  H = ( Proj_W[ Z F^{23} ], Proj_W[ Z F^{31} ], Proj_W[ Z F^{12} ] )")
print("These need not coincide with the measured fields (E,B).")

# Antisymmetric source-coupled tensor components
G10, G20, G30 = D1, D2, D3
G01, G02, G03 = -D1, -D2, -D3
G23, G31, G12 = H1, H2, H3
G32, G13, G21 = -H1, -H2, -H3

# nu = 0 (Gauss-like law)
lhs0 = sp.diff(G10, x) + sp.diff(G20, y) + sp.diff(G30, z) + L0 + G0
rhs0 = mu0 * rho
print("\nGauss-like law from the projected inhomogeneous equation (nu = 0):")
print("  ", sp.sstr(sp.Eq(lhs0, rhs0)))
print("  i.e. div D + Leak0 + Gauge0 = mu0 rho_proj")

# nu = 1,2,3 (Ampere-like law)
lhs1 = sp.diff(G01, t) + sp.diff(G21, y) + sp.diff(G31, z) + L1 + G1
lhs2 = sp.diff(G02, t) + sp.diff(G32, z) + sp.diff(G12, x) + L2 + G2
lhs3 = sp.diff(G03, t) + sp.diff(G13, x) + sp.diff(G23, y) + L3 + G3

curlH1 = sp.diff(H3, y) * (-1) + sp.diff(H2, z)
curlH2 = sp.diff(H1, z) + sp.diff(H3, x) * (-1)
curlH3 = sp.diff(H2, x) + sp.diff(H1, y) * (-1)

amp1_target = (sp.diff(H2, z) - sp.diff(H3, y)) - sp.diff(D1, t) + L1 + G1
amp2_target = (sp.diff(H3, x) - sp.diff(H1, z)) - sp.diff(D2, t) + L2 + G2
amp3_target = (sp.diff(H1, y) - sp.diff(H2, x)) - sp.diff(D3, t) + L3 + G3

check_amp = [sp.simplify(lhs1 - amp1_target), sp.simplify(lhs2 - amp2_target), sp.simplify(lhs3 - amp3_target)]

print("\nAmpere-like components from the projected inhomogeneous equation:")
print("  ", sp.sstr(sp.Eq(amp1_target, mu0 * J1)))
print("  ", sp.sstr(sp.Eq(amp2_target, mu0 * J2)))
print("  ", sp.sstr(sp.Eq(amp3_target, mu0 * J3)))
print("\nVerification residues (all should be 0):")
print("  ", check_amp)
for i, residue in enumerate(check_amp, start=1):
    assert_zero(f"Ampere component {i}", residue)

print("\nCompact inhomogeneous vector form:")
print("  div D + Leak0 + Gauge0 = mu0 rho_proj")
print("  curl H - partial_t D + Leak_vec + Gauge_vec = mu0 J_proj")


# -----------------------------------------------------------------------------
# Physical interpretation block
# -----------------------------------------------------------------------------
line("3) Physical interpretation")
print("Measured fields are built from Proj_W[F_{MN}] and satisfy the homogeneous equations directly.")
print("Source-coupled flux fields are built from Proj_W[Z F^{MN}] and enter the inhomogeneous equations.")
print("Therefore a projected 4+1 Maxwell theory naturally has two field layers:")
print("  (E_meas, B_meas)     -> what the brane observer reports")
print("  (D_flux, H_flux)     -> what the projected inhomogeneous operator couples to")
print()
print("They coincide only under additional closure assumptions, such as trivial Z(w) structure or a dedicated reduction ansatz.")
print("The extra terms Leak_mu are generated by Z F^{w mu} and encode exchange with the hidden w-direction.")
print("The Gauge_mu terms are the projected gauge-driver contributions.")


# -----------------------------------------------------------------------------
# Compact summary
# -----------------------------------------------------------------------------
line("4) Compact summary")
print("Homogeneous projected system:")
print("  div B_meas = 0")
print("  curl E_meas + partial_t B_meas = 0")
print()
print("Inhomogeneous projected system:")
print("  div D_flux + Leak0 + Gauge0 = mu0 rho_proj")
print("  curl H_flux - partial_t D_flux + Leak_vec + Gauge_vec = mu0 J_proj")
print()
print("So the projected theory looks like an open-system effective medium, not a closed copy of reduced 3+1 Maxwell.")


# -----------------------------------------------------------------------------
# Concrete projection audit
# -----------------------------------------------------------------------------
line("5) Concrete projection audit")

w = sp.symbols("w", real=True)
Wg = sp.exp(-w**2) / sp.sqrt(sp.pi)
Wgp = sp.diff(Wg, w)
Zg = sp.exp(-w**2)

def Pg(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.integrate(Wg * expr, (w, -sp.oo, sp.oo)))


def Pgp(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.integrate(Wgp * expr, (w, -sp.oo, sp.oo)))


def boundary_value(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.limit(expr, w, sp.oo) - sp.limit(expr, w, -sp.oo))


F10_bulk = 1 + w**2
E_meas_example = Pg(F10_bulk)
D_flux_example = Pg(Zg * F10_bulk)
assert_nonzero("measured and flux fields differ for nontrivial Z", E_meas_example - D_flux_example)

# Explicit bulk potential to make the projected Bianchi and inhomogeneous
# vector laws concrete instead of purely notational.
A0_bulk = x * z * (1 + w**2)
A1_bulk = t * y * (1 + w**2)
A2_bulk = t * z * (1 + w**2)
A3_bulk = x * y * (1 + w**2)

F10_bulk = sp.diff(A0_bulk, x) - sp.diff(A1_bulk, t)
F20_bulk = sp.diff(A0_bulk, y) - sp.diff(A2_bulk, t)
F30_bulk = sp.diff(A0_bulk, z) - sp.diff(A3_bulk, t)
F23_bulk = sp.diff(A3_bulk, y) - sp.diff(A2_bulk, z)
F31_bulk = sp.diff(A1_bulk, z) - sp.diff(A3_bulk, x)
F12_bulk = sp.diff(A2_bulk, x) - sp.diff(A1_bulk, y)

E1_bulk_proj, E2_bulk_proj, E3_bulk_proj = [Pg(expr) for expr in (F10_bulk, F20_bulk, F30_bulk)]
B1_bulk_proj, B2_bulk_proj, B3_bulk_proj = [Pg(expr) for expr in (F23_bulk, F31_bulk, F12_bulk)]
assert_zero("concrete projected div B", sp.diff(B1_bulk_proj, x) + sp.diff(B2_bulk_proj, y) + sp.diff(B3_bulk_proj, z))
assert_zero(
    "concrete projected Faraday 1",
    sp.diff(B1_bulk_proj, t) + sp.diff(E3_bulk_proj, y) - sp.diff(E2_bulk_proj, z),
)
assert_zero(
    "concrete projected Faraday 2",
    sp.diff(B2_bulk_proj, t) + sp.diff(E1_bulk_proj, z) - sp.diff(E3_bulk_proj, x),
)
assert_zero(
    "concrete projected Faraday 3",
    sp.diff(B3_bulk_proj, t) + sp.diff(E2_bulk_proj, x) - sp.diff(E1_bulk_proj, y),
)

Fw1 = w
ZFw1 = Zg * Fw1
leak1 = -Pgp(ZFw1)
projected_transverse_derivative = Pg(sp.diff(ZFw1, w))
boundary_ZFw1 = boundary_value(Wg * ZFw1)
assert_zero("localized leakage boundary term", boundary_ZFw1)
assert_zero("leakage IBP with explicit boundary", projected_transverse_derivative - (boundary_ZFw1 - Pgp(ZFw1)))
assert_zero("leakage equals projected transverse derivative", leak1 - projected_transverse_derivative)
assert_zero("nonzero vector leakage value", leak1 - sp.sqrt(2) / 4)
assert_nonzero("mutated leakage IBP sign should fail", projected_transverse_derivative - (boundary_ZFw1 + Pgp(ZFw1)))

Fw0_bulk = sp.diff(A0_bulk, w)
Fw1_bulk = sp.diff(A1_bulk, w)
Fw2_bulk = sp.diff(A2_bulk, w)
Fw3_bulk = sp.diff(A3_bulk, w)
D1_bulk_proj, D2_bulk_proj, D3_bulk_proj = [Pg(Zg * expr) for expr in (F10_bulk, F20_bulk, F30_bulk)]
H1_bulk_proj, H2_bulk_proj, H3_bulk_proj = [Pg(Zg * expr) for expr in (F23_bulk, F31_bulk, F12_bulk)]
boundary_flux_terms = [boundary_value(Wg * Zg * expr) for expr in (Fw0_bulk, Fw1_bulk, Fw2_bulk, Fw3_bulk)]
for idx, value in enumerate(boundary_flux_terms):
    assert_zero(f"concrete transverse boundary term {idx}", value)
leak0_bulk, leak1_bulk, leak2_bulk, leak3_bulk = [-Pgp(Zg * expr) for expr in (Fw0_bulk, Fw1_bulk, Fw2_bulk, Fw3_bulk)]
for idx, expr in enumerate((Fw0_bulk, Fw1_bulk, Fw2_bulk, Fw3_bulk)):
    assert_zero(
        f"concrete transverse IBP identity {idx}",
        Pg(sp.diff(Zg * expr, w)) - (boundary_flux_terms[idx] - Pgp(Zg * expr)),
    )
J0_bulk = (
    sp.diff(Zg * F10_bulk, x) + sp.diff(Zg * F20_bulk, y) + sp.diff(Zg * F30_bulk, z) + sp.diff(Zg * Fw0_bulk, w)
) / mu0
J1_bulk = (
    sp.diff(-Zg * F10_bulk, t) + sp.diff(-Zg * F12_bulk, y) + sp.diff(Zg * F31_bulk, z) + sp.diff(Zg * Fw1_bulk, w)
) / mu0
J2_bulk = (
    sp.diff(-Zg * F20_bulk, t) + sp.diff(-Zg * F23_bulk, z) + sp.diff(Zg * F12_bulk, x) + sp.diff(Zg * Fw2_bulk, w)
) / mu0
J3_bulk = (
    sp.diff(-Zg * F30_bulk, t) + sp.diff(-Zg * F31_bulk, x) + sp.diff(Zg * F23_bulk, y) + sp.diff(Zg * Fw3_bulk, w)
) / mu0
assert_zero(
    "concrete projected Gauss law",
    sp.diff(D1_bulk_proj, x) + sp.diff(D2_bulk_proj, y) + sp.diff(D3_bulk_proj, z) + leak0_bulk - mu0 * Pg(J0_bulk),
)
assert_zero(
    "concrete projected Ampere 1",
    (sp.diff(H2_bulk_proj, z) - sp.diff(H3_bulk_proj, y)) - sp.diff(D1_bulk_proj, t) + leak1_bulk - mu0 * Pg(J1_bulk),
)
assert_zero(
    "concrete projected Ampere 2",
    (sp.diff(H3_bulk_proj, x) - sp.diff(H1_bulk_proj, z)) - sp.diff(D2_bulk_proj, t) + leak2_bulk - mu0 * Pg(J2_bulk),
)
assert_zero(
    "concrete projected Ampere 3",
    (sp.diff(H1_bulk_proj, y) - sp.diff(H2_bulk_proj, x)) - sp.diff(D3_bulk_proj, t) + leak3_bulk - mu0 * Pg(J3_bulk),
)
assert_nonzero(
    "mutated concrete Faraday sign should fail",
    sp.diff(B1_bulk_proj, t) - sp.diff(E3_bulk_proj, y) + sp.diff(E2_bulk_proj, z),
)

print("Concrete checks: Faraday/Ampere sign maps, explicit boundary discharge, E-vs-D distinction, and nonzero leakage all pass.")
print("STATUS: PASS")
