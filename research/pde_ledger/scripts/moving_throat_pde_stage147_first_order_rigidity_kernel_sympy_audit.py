#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t):
    print("\n" + "="*88)
    print(t)
    print("="*88)

Pi = sp.symbols("Pi", positive=True, real=True)
x, eps = sp.symbols("x eps", real=True)
kap = sp.pi/2
c = sp.cos(sp.pi*x/2)
Kq = sp.cosh(kap*(1-x))/sp.cosh(kap)

gPi = 2*Pi*(2*Pi*sp.exp(Pi)+sp.pi)/((4*Pi**2+sp.pi**2)*(sp.exp(Pi)-1))
Sformula = sp.simplify(
    Pi*(kap*sp.tanh(kap) + Pi*(sp.exp(-Pi)*sp.sech(kap) - 1))
    / ((1-sp.exp(-Pi))*(kap**2-Pi**2))
)

rF1 = sp.sqrt(12*sp.Rational(37,20)**2/sp.pi**2 - 1)
gminus = sp.simplify(rF1 - sp.sqrt(1+rF1**2)/2)
Pi_star = sp.N(sp.nsolve(gPi - gminus, 1.5), 40)
g_star = sp.N(gPi.subs(Pi, Pi_star), 40)
S_star = sp.N(Sformula.subs(Pi, Pi_star), 40)
gp_star = sp.N(sp.diff(gPi, Pi).subs(Pi, Pi_star), 40)
Sp_star = sp.N(sp.diff(Sformula, Pi).subs(Pi, Pi_star), 40)

Sigma_star = sp.N(Pi_star/(1-S_star/4), 40)
T_star = sp.N(sp.sqrt(9*Sigma_star/20), 40)

AT = sp.N(
    -(sp.Rational(9,1)/(40*T_star)) * (
        1/(gp_star*(1-S_star/4)) + Pi_star*Sp_star/(4*gp_star*(1-S_star/4)**2)
    ),
    30
)
BT = sp.N(
    (sp.Rational(9,1)/(40*T_star)) * Pi_star/(4*(1-S_star/4)**2),
    30
)

banner("Stage 147 audit: first-order rigidity coefficients")
print("Pi_* =", Pi_star)
print("Sigma_* =", Sigma_star)
print("T_* =", T_star)
print("A_T =", AT)
print("B_T =", BT)
print("|A_T|/B_T =", sp.N(abs(AT)/BT, 20))

# Weight kernel representation
gbar, Sbar = sp.symbols("gbar Sbar", real=True)
dT = sp.simplify(eps*(AT*(gbar-gminus) + BT*(Sbar-Sformula.subs(Pi, Pi_star))))
print("delta T_m =")
sp.pprint(dT)

# Centered kernel form: constant pieces cancel against normalization.
Wcenter = sp.simplify(AT*(c-gminus) + BT*(Kq-Sformula.subs(Pi, Pi_star)))
print("Centered rigidity kernel W_*(x) =")
sp.pprint(Wcenter)

print("\nStage 147 complete.")
