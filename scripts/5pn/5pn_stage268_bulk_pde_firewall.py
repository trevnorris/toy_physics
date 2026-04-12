#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

"""
Stage 268 — bulk-PDE firewall check.

Purpose
-------
Formalize the reduced statement that, once the lower-order hierarchy is frozen,
the remaining 5PN / 2.5PN / 4PN endgame lives in moving-throat branch data,
not in a retuning of the parent GNLS bulk medium.

This is not a proof of the full moving-throat PDE realization theorem.
It is the exact reduced bookkeeping statement we should now enforce.
"""

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

banner("STAGE 268 — BULK-PDE FIREWALL CHECK")

subbanner("I. Lower-order hierarchy stays frozen")
n = sp.symbols("n", positive=True, real=True)
kappa_rho, kappa_add, kappa_PV, beta_1PN = sp.symbols(
    "kappa_rho kappa_add kappa_PV beta_1PN", real=True
)

frozen = {
    n: sp.Integer(5),
    kappa_rho: sp.Integer(1),
    kappa_add: sp.Rational(1, 2),
    kappa_PV: sp.Rational(3, 2),
    beta_1PN: sp.Integer(3),
}
print("Frozen carry-forward hierarchy:")
for k, v in frozen.items():
    print(f"  {sp.sstr(k)} = {sp.sstr(v)}")

subbanner("II. Conservative pole-fraction firewall")
eps2, eps4, chi, a2, a4 = sp.symbols("eps_2 eps_4 chi a_2 a_4", real=True)
c_pole = sp.simplify((1 + eps4) / (4 * (1 + eps2) ** 2))
print("Exact grouped-P2 conservative pole fraction:")
sp.pprint(sp.Eq(sp.Symbol("c_pole"), c_pole))

c_pole_iso = sp.simplify(c_pole.subs({eps2: 0, eps4: 0}))
print("Isotropic branch value:")
sp.pprint(sp.Eq(sp.Symbol("c_pole"), c_pole_iso))

c_pole_weak = sp.series(
    c_pole.subs({eps2: a2 * chi**2, eps4: a4 * chi**2}),
    chi,
    0,
    4,
).removeO()
print("First weak-mixing correction with eps_2=a_2 chi^2, eps_4=a_4 chi^2:")
sp.pprint(sp.Eq(sp.Symbol("c_pole"), sp.simplify(c_pole_weak)))

subbanner("III. Remaining reduced closure packet")
a2b, b2, a4b, b4, aP0, bP0 = sp.symbols("a_2 b_2 a_4 b_4 a_P0 b_P0", real=True)
Delta_pole, Delta_norm = sp.symbols("Delta_pole Delta_norm", real=True)
qtr, qnt, qeta = sp.symbols("q_tr q_nt q_eta", real=True)

Delta_branch = sp.Matrix([a2b, b2, a4b, b4, aP0, bP0, Delta_pole, Delta_norm])
Delta_orbit = sp.Matrix([qtr, qnt, qeta])

print("Reduced closure data packet:")
print("Delta_branch =")
sp.pprint(Delta_branch)
print("Delta_orbit =")
sp.pprint(Delta_orbit)

subbanner("IV. Firewall interpretation")
print("Once the lower-order hierarchy is frozen, the reduced endgame depends on:")
print("  (1) grouped branch packet Delta_branch")
print("  (2) orbit/selector packet Delta_orbit")
print("and not on a new retuning of the parent EOS exponent n.")
print()
print("Working rule:")
print("  keep the parent GNLS + localized Maxwell bulk theory fixed;")
print("  continue only on the moving-throat branch / boundary-response side.")
