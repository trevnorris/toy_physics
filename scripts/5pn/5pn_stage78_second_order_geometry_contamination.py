
#!/usr/bin/env python3
"""
5pn_stage78_second_order_geometry_contamination.py

Stage 78 audit: the first nonzero geometry contamination appears only at second
order in anisotropy/mixing.
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

banner("STAGE 78 — FIRST NONZERO GEOMETRY CONTAMINATION IS SECOND ORDER")

omega = sp.symbols("omega", real=True)
chi, M0 = sp.symbols("chi M0", real=True)
Kstat, Kpole, Omega_Q = sp.symbols("K_stat K_pole Omega_Q", positive=True, real=True)
G0, G2, G4 = sp.symbols("G0 G2 G4", positive=True, real=True)
e2hat, e4hat = sp.symbols("e2hat e4hat", real=True)

subbanner("78.1 — Minimal mixed scalar-geometry / grouped-P2 model")
D_q = sp.simplify(Kstat + Kpole / (1 - omega**2 / Omega_Q**2))
D_g = sp.simplify(G0 + G2 * omega**2 + G4 * omega**4)
D_eff = sp.simplify(D_q - chi**2 * M0**2 / D_g)
print("D_eff(omega) =")
sp.pprint(D_eff)

subbanner("78.2 — Exact low-frequency contamination coefficients")
inv_Dg = sp.series(1 / D_g, omega, 0, 6).removeO()
print("1/D_g(omega) =")
sp.pprint(sp.expand(inv_Dg))

K_g2_eff = sp.simplify(chi**2 * M0**2 * G2 / G0**2)
K_g4_eff = sp.simplify(chi**2 * M0**2 * (G0 * G4 - G2**2) / G0**3)
print("K_(g,2)^eff =", K_g2_eff)
print("K_(g,4)^eff =", K_g4_eff)

eps2 = sp.simplify(Omega_Q**2 * K_g2_eff / Kpole)
eps4 = sp.simplify(Omega_Q**4 * K_g4_eff / Kpole)
print("eps_2 =", eps2)
print("eps_4 =", eps4)

subbanner("78.3 — Contact/pole fraction deviation is O(chi^2)")
c_pole_abstract = sp.simplify((1 + chi**2 * e4hat) / (4 * (1 + chi**2 * e2hat)**2))
series_c2 = sp.series(c_pole_abstract, chi, 0, 3).removeO()
expected = sp.simplify(sp.Rational(1, 4) * (1 + chi**2 * (e4hat - 2 * e2hat)))
expect_zero("c_pole expansion", sp.expand(series_c2 - expected))

banner("STAGE 78 FINAL LEDGER")
print("The first nonzero geometry contamination requires an explicit l=0 <-> l=2 mixing")
print("source and appears only at O(chi^2):")
print("  eps_2, eps_4 = O(chi^2),")
print("so the deviation from the Stage-74 3/4 + 1/4 split also begins only at O(chi^2).")
