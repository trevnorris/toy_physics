#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t):
    print("\n" + "="*88)
    print(t)
    print("="*88)

Pi = sp.symbols("Pi", positive=True, real=True)
x, eps, lam = sp.symbols("x eps lam", real=True)
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
# Paper-literal anchor (EXTERNAL, appendix part04:846/848): AT/BT reproduce the published
# A_T/B_T. A wrong AT formula (e.g. a dropped chain term) FAILS here.
assert abs(sp.N(AT, 30) - sp.Float("-4.27263956256927")) < sp.Float("1e-11"), f"AT vs paper A_T: {AT}"
assert abs(sp.N(BT, 30) - sp.Float("0.134875005736706")) < sp.Float("1e-11"), f"BT vs paper B_T: {BT}"

banner("Stage 148 audit: representative non-exponential families")

# Uniform
g_u = sp.N(2/sp.pi, 30)
S_u = sp.N(2*sp.tanh(sp.pi/2)/sp.pi, 30)

dPi_u = sp.N(-(g_u-g_star)/gp_star, 30)
dT_u = sp.N(AT*(g_u-g_star) + BT*(S_u-S_star), 30)

print("uniform: g_u =", g_u)
print("uniform: S_u =", S_u)
print("uniform: dPi/eps =", dPi_u)
print("uniform: dT/eps =", dT_u)

# Self-matched derivative family
g_d = sp.N(sp.pi/4, 30)
S_d = sp.N((1+sp.sinh(sp.pi/2))/(2*sp.cosh(sp.pi/2)), 30)

dPi_d = sp.N(-(g_d-g_star)/gp_star, 30)
dT_d = sp.N(AT*(g_d-g_star) + BT*(S_d-S_star), 30)

print("derivative: g_d =", g_d)
print("derivative: S_d =", S_d)
print("derivative: dPi/eps =", dPi_d)
print("derivative: dT/eps =", dT_d)

# Convex interpolation
g_lam = sp.simplify((1-lam)*g_u + lam*g_d)
S_lam = sp.simplify((1-lam)*S_u + lam*S_d)
dPi_lam = sp.simplify(-(g_lam-g_star)/gp_star)
dT_lam = sp.simplify(AT*(g_lam-g_star) + BT*(S_lam-S_star))

print("dPi_lambda/eps =")
sp.pprint(sp.expand(dPi_lam))
print("dT_lambda/eps =")
sp.pprint(sp.expand(dT_lam))

lam_Pi_zero = sp.N(sp.solve(sp.Eq(dPi_lam,0), lam)[0], 30)
lam_T_zero = sp.N(sp.solve(sp.Eq(dT_lam,0), lam)[0], 30)
print("lambda_(Pi,0) =", lam_Pi_zero)
print("lambda_(T,0) =", lam_T_zero)
print("1 - lambda_(Pi,0) =", sp.N(1-lam_Pi_zero, 30))

# Stage 126 positive-family compensation closed form (see notes section 3).
# GUARD: 100 under the radical is FORCED by rF1 (12*(37/20)^2 = 4107/100, so
# rF1^2 = (4107 - 100*pi^2)/(100*pi^2)); do NOT substitute any other constant.
xi_star_closed = (-37*sp.sqrt(3) - 5*sp.pi**2 + 2*sp.sqrt(4107 - 100*sp.pi**2)) / (5*(8 - sp.pi**2))
xi_star = sp.N(xi_star_closed, 30)
print("xi_* (Stage 126 closed form) =", xi_star)
# Preferred (EXACT, no Pi_star / no nsolve): 1 - lambda_(Pi,0) has the closed form
# (pi/4 - gminus)/(pi/4 - 2/pi), since the convex family is g_lam = (1-lam)*(2/pi)+lam*(pi/4)
# and lambda_(Pi,0) solves g_lam = g_* = gminus. Build it from the EXACT gminus (NOT the
# sp.N-numericized g_star/lam_Pi_zero, which would defeat the symbolic reduction). Note
# 1 + rF1^2 = 4107/(100 pi^2) so sqrt(1+rF1^2) = 37*sqrt(3)/(10 pi), and sqrt(4107)=37*sqrt(3),
# so SymPy collapses the nested radicals and reduces the difference to exact 0.
rF1_exact = sp.sqrt(12*sp.Rational(37, 20)**2/sp.pi**2 - 1)
gminus_exact = rF1_exact - sp.sqrt(1 + rF1_exact**2)/2
one_minus_lam_exact = (sp.pi/4 - gminus_exact) / (sp.pi/4 - 2/sp.pi)
exact_resid = sp.simplify(one_minus_lam_exact - xi_star_closed)
print("exact (1-lambda_(Pi,0)) - xi_* =", exact_resid)
residual = sp.N((1 - lam_Pi_zero) - xi_star, 30)   # numeric, for transcript display
print("(1-lambda_(Pi,0)) - xi_* =", residual)
assert exact_resid == 0, f"Stage 148 D4 consistency (exact) failed: residual = {exact_resid}"

print("\nStage 148 complete.")
