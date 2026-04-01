#!/usr/bin/env python3
"""
moving_throat_pde_stage53_gnls_wall_shell_sympy_audit.py

SymPy audit for Stage 53:
explicit GNLS wall-shell reduction of the first support/source branch.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 53 — EXPLICIT GNLS WALL-SHELL REDUCTION")

a, L, ell = sp.symbols("a L ell", positive=True, real=True)
rho_w, c_sw, V0 = sp.symbols("rho_w c_sw V0", positive=True, real=True)
m, hbar = sp.symbols("m hbar", positive=True, real=True)
If, Ig = sp.symbols("I_f I_g", positive=True, real=True)
pi = sp.pi

Hw = sp.simplify(m * c_sw**2 / rho_w)

Nphiphi = sp.simplify(4 * pi * a**2 * ell * If)
Gphiphi = sp.simplify(4 * pi * a**2 * Ig / ell)

Tx = sp.simplify(hbar**2 * Nphiphi / (4 * m * rho_w))
Kx = sp.simplify(Hw * Nphiphi + hbar**2 * Gphiphi / (4 * m * rho_w))
kappa = sp.simplify(Kx * L**2 / Tx)

print("N_(phi phi) =", Nphiphi)
print("G_(phi phi) =", Gphiphi)
print("T_X =", Tx)
print("K_X =", Kx)
print("kappa =", kappa)

kappa_expected = sp.simplify(4 * (m * c_sw * L / hbar)**2 + (Ig / If) * (L / ell)**2)
expect_zero("kappa - expected", kappa - kappa_expected)

J1 = sp.simplify(If / Hw)
Wwall = sp.simplify(4 * pi * a**2 * L**2 * J1 * V0**2 / (Tx * ell))
Wwall_expected = sp.simplify(4 * rho_w**2 * V0**2 * L**2 / (hbar**2 * c_sw**2 * ell**2))
print("J_1 =", J1)
print("W_wall =", Wwall)
expect_zero("W_wall - expected", Wwall - Wwall_expected)

gphi = sp.simplify(V0 / ell)
I1 = sp.simplify(Nphiphi / Hw)
Xi = sp.simplify(gphi**2 * I1 * L**2 / Tx)
print("g_phi =", gphi)
print("Xi =", Xi)
expect_zero("Xi - W_wall", Xi - Wwall)

banner("STAGE 53 THEOREM LEDGER")
print("T_X = pi a^2 ell I_f hbar^2 / (m rho_w)")
print("K_X = 4 pi a^2 ell I_f H_w + pi a^2 I_g hbar^2 / (m rho_w ell)")
print("kappa = 4 (m c_(s,w) L / hbar)^2 + (I_g/I_f)(L/ell)^2")
print("W_wall = Xi = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2)")
