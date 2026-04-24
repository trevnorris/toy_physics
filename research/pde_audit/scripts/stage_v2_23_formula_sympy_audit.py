#!/usr/bin/env python3
"""Stage V2-23 formula-only SymPy audit."""
from __future__ import annotations

import sympy as sp

w = sp.symbols("w")
D0, D2, D4, N0, N2, N4 = sp.symbols("D0 D2 D4 N0 N2 N4", nonzero=True)
B4, Z4, M, B2, Z2 = sp.symbols("B4 Z4 M B2 Z2", nonzero=True)
a, cs, G, c, mhat, Sport = sp.symbols("a cs G c mhat Sport", positive=True)

D = D0 + D2*w**2 + D4*w**4
N = N0 + N2*w**2 + N4*w**4

Y = sp.series(D0/D, w, 0, 6).removeO().expand()
u2 = -D2/D0
u4 = (D2**2 - D0*D4)/D0**2

Pref = sp.series(D0*N/D**2, w, 0, 6).removeO().expand()
P0 = N0/D0
P2 = (D0*N2 - 2*D2*N0)/D0**2
P4 = (D0**2*N4 - 2*D0*(D2*N2 + D4*N0) + 3*D2**2*N0)/D0**3

# One-pole condition in the isotropic bundle:
# D2 = -(M+B2+Z2), D4 = -(B4+Z4).
A = M + B2 + Z2
C = B4 + Z4
R_pole = D0*C - 3*A**2
one_pole_equiv = sp.simplify((u4 - 4*u2**2).subs({D2: -A, D4: -C}) * D0**2 - R_pole)

N2_const = 2*D2*N0/D0
N4_const = N0*(D2**2 + 2*D0*D4)/D0**2

P0target = 54*G*cs**5/(5*a**5*c**5)
gamma_eff = mhat**2*Sport*P0target*a**5/(27*cs**5)
gamma_GR = 2*G/(5*c**5)

checks = {
    "response_u2": sp.simplify(Y.coeff(w, 2) - u2) == 0,
    "response_u4": sp.simplify(Y.coeff(w, 4) - u4) == 0,
    "prefactor_P0": sp.simplify(Pref.coeff(w, 0) - P0) == 0,
    "prefactor_P2": sp.simplify(Pref.coeff(w, 2) - P2) == 0,
    "prefactor_P4": sp.simplify(Pref.coeff(w, 4) - P4) == 0,
    "one_pole_equivalence": one_pole_equiv == 0,
    "constant_prefactor_P2": sp.simplify(P2.subs(N2, N2_const)) == 0,
    "constant_prefactor_P4": sp.simplify(P4.subs({N2: N2_const, N4: N4_const})) == 0,
    "normalization_equivalence": sp.simplify(gamma_eff.subs({mhat: 1, Sport: 1}) - gamma_GR) == 0,
}

lines = ["Stage V2-23 formula-only SymPy audit", "=" * 48]
for name, ok in checks.items():
    lines.append(f"{name}: {'PASS' if ok else 'FAIL'}")
lines.append(f"checks_passed: {sum(bool(v) for v in checks.values())}")
lines.append(f"checks_total: {len(checks)}")

output = "\n".join(lines) + "\n"
print(output)
