#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t):
    print("\n" + "="*88)
    print(t)
    print("="*88)

def expect_zero(name, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

Pi = sp.symbols("Pi", positive=True, real=True)
x = sp.symbols("x", real=True)
kap = sp.pi/2

# Exact exponential-source overlap and mixed-response functions.
Sigma = Pi*sp.exp(-Pi*x)/(1-sp.exp(-Pi))
Kq = sp.cosh(kap*(1-x))/sp.cosh(kap)

gPi = 2*Pi*(2*Pi*sp.exp(Pi)+sp.pi)/((4*Pi**2+sp.pi**2)*(sp.exp(Pi)-1))
Sformula = sp.simplify(
    Pi*(kap*sp.tanh(kap) + Pi*(sp.exp(-Pi)*sp.sech(kap) - 1))
    / ((1-sp.exp(-Pi))*(kap**2-Pi**2))
)

banner("Stage 146 audit: exact source moments")
print("g(Pi) =")
sp.pprint(gPi)
print("S_q(Pi) =")
sp.pprint(Sformula)

# Numerical kernel cross-check against direct integration.
for val in [sp.Integer(1), sp.Rational(3,2), sp.Rational(5,2)]:
    num_int = sp.N(sp.integrate((Sigma*Kq).subs(Pi, val), (x,0,1)), 30)
    num_formula = sp.N(Sformula.subs(Pi, val), 30)
    diff = sp.N(num_int - num_formula, 20)
    print(f"kernel check at Pi={val}: integral={num_int}, formula={num_formula}, diff={diff}")
    if abs(float(diff)) > 1e-12:
        raise AssertionError("Kernel formula mismatch.")

# Family-1 canonical point
rF1 = sp.sqrt(12*sp.Rational(37,20)**2/sp.pi**2 - 1)
gminus = sp.simplify(rF1 - sp.sqrt(1+rF1**2)/2)
Pi_star = sp.N(sp.nsolve(gPi - gminus, 1.5), 30)
print("Pi_* =", Pi_star)

g_star = sp.N(gPi.subs(Pi, Pi_star), 30)
S_star = sp.N(Sformula.subs(Pi, Pi_star), 30)
gp_star = sp.N(sp.diff(gPi, Pi).subs(Pi, Pi_star), 30)
Sp_star = sp.N(sp.diff(Sformula, Pi).subs(Pi, Pi_star), 30)
print("g_* =", g_star)
print("S_* =", S_star)
print("g_*' =", gp_star)
print("S_*' =", Sp_star)

# Generic positive convex deformation moments
gbar, Sbar, eps = sp.symbols("gbar Sbar eps", real=True)
dPi = sp.simplify(-eps*(gbar-gminus)/sp.diff(gPi, Pi).subs(Pi, Pi_star))
dS = sp.simplify(eps*(Sbar - Sformula.subs(Pi, Pi_star)) + sp.diff(Sformula, Pi).subs(Pi, Pi_star)*dPi)

banner("First-order canonical retuning law")
print("delta Pi =")
sp.pprint(dPi)
print("delta S =")
sp.pprint(sp.simplify(dS))

# Exact affine law for convex family moments.
g_eps = (1-eps)*gminus + eps*gbar
S_eps = (1-eps)*Sformula.subs(Pi, Pi_star) + eps*Sbar
expect_zero("g_eps affine law", sp.expand(g_eps - (gminus + eps*(gbar-gminus))))
expect_zero("S_eps affine law", sp.expand(S_eps - (Sformula.subs(Pi, Pi_star) + eps*(Sbar-Sformula.subs(Pi, Pi_star)))))

print("\nStage 146 complete.")
