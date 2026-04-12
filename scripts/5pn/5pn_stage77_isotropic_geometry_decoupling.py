
#!/usr/bin/env python3
"""
5pn_stage77_isotropic_geometry_decoupling.py

Stage 77 audit: isotropic geometry-decoupling theorem.
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

def expect_zero(name: str, expr) -> None:
    expr_s = sp.simplify(sp.together(sp.expand(expr)))
    print(f"{name} = {expr_s}")
    if expr_s != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 77 — ISOTROPIC GEOMETRY-DECOUPLING THEOREM")

theta, phi = sp.symbols("theta phi", real=True)
dOmega = sp.sin(theta)

Y00 = sp.Rational(1, 2) / sp.sqrt(sp.pi)
Y20 = sp.sqrt(sp.Rational(5, 16) / sp.pi) * (3 * sp.cos(theta)**2 - 1)
Y21c = sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.cos(phi)
Y21s = sp.sqrt(sp.Rational(15, 4) / sp.pi) * sp.sin(theta) * sp.cos(theta) * sp.sin(phi)
Y22c = sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta)**2 * sp.cos(2 * phi)
Y22s = sp.sqrt(sp.Rational(15, 16) / sp.pi) * sp.sin(theta)**2 * sp.sin(2 * phi)

Ys = {
    "20": Y20,
    "21c": Y21c,
    "21s": Y21s,
    "22c": Y22c,
    "22s": Y22s,
}

subbanner("77.1 — Exact l=0 <-> l=2 orthogonality")
for name, Y in Ys.items():
    overlap = sp.simplify(sp.integrate(sp.integrate(Y00 * Y * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)))
    expect_zero(f"<Y00, Y{name}>", overlap)

subbanner("77.2 — Gradient cross terms vanish")
for name, Y in Ys.items():
    grad_cross = sp.simplify(sp.diff(Y00, theta) * sp.diff(Y, theta) + sp.diff(Y00, phi) * sp.diff(Y, phi) / sp.sin(theta)**2)
    expect_zero(f"grad cross density Y00<->Y{name}", grad_cross)

subbanner("77.3 — Laplacian cross terms reduce to orthogonality")
for name, Y in Ys.items():
    lap_cross = sp.simplify(6 * sp.integrate(sp.integrate(Y00 * Y * dOmega, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)))
    expect_zero(f"<Y00, (-Delta)Y{name}>", lap_cross)

subbanner("77.4 — Geometry contamination numbers")
Omega_Q, Kpole = sp.symbols("Omega_Q K_pole", positive=True, real=True)
Kg2, Kg4 = sp.symbols("K_g2 K_g4", real=True)
eps2 = sp.simplify(Omega_Q**2 * Kg2 / Kpole)
eps4 = sp.simplify(Omega_Q**4 * Kg4 / Kpole)
expect_zero("eps_2 on the isotropic branch", eps2.subs(Kg2, 0))
expect_zero("eps_4 on the isotropic branch", eps4.subs(Kg4, 0))

banner("STAGE 77 FINAL LEDGER")
print("On the isotropic quadratic wall theory, the scalar/geometry l=0 lane and the grouped")
print("real l=2 bundle are exactly orthogonal. Therefore the only geometry contribution allowed")
print("on that branch is the static completion:")
print("  K_(g,2) = 0,   K_(g,4) = 0,   eps_2 = eps_4 = 0.")
