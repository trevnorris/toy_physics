#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t):
    line="="*88
    print("\n"+line); print(t); print(line)

def subbanner(t):
    line="-"*88
    print("\n"+line); print(t); print(line)

def expect_zero(name, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

Pi = sp.symbols("Pi", positive=True, real=True)
pi = sp.pi

banner("STAGE 125 — SELF-CONSISTENT MOUTH-BRANCH LAW")

r = sp.sqrt(sp.Integer(4107) - 100*pi**2)/(10*pi)
gPi = 2*Pi*(2*Pi*sp.exp(Pi)+pi)/((4*Pi**2+pi**2)*(sp.exp(Pi)-1))
Sq = Pi*(((pi/2)*sp.tanh(pi/2)) + Pi*(sp.exp(-Pi)*sp.sech(pi/2)-1))/((1-sp.exp(-Pi))*(((pi/2)**2)-Pi**2))
Rq = sp.simplify((gPi-r)**2/(1+r**2))
Sigma0 = Pi/(1-Rq*Sq)
That = sp.sqrt(sp.Rational(9,20)*Sigma0)

subbanner("Core-to-mouth reduction")
print("r_F1 =", r)
print("g_Pi =", gPi)
print("S_q(Pi) =", Sq)
print("R_q(Pi) =", Rq)
print("Sigma0(Pi) =", Sigma0)
print("That(Pi) =", That)

gminus = sp.simplify(r - sp.Rational(1,2)*sp.sqrt(1+r**2))
Rq_minus = sp.simplify(((gminus-r)**2)/(1+r**2))
expect_zero("R_q(g_minus)-1/4", Rq_minus - sp.Rational(1,4))

subbanner("Canonical compensation point")
Pi_star = sp.N(sp.nsolve(sp.Eq(gPi, gminus), 1.5), 30)
g_star = sp.N(gPi.subs(Pi, Pi_star), 30)
Rq_star = sp.N(Rq.subs(Pi, Pi_star), 30)
Sq_star = sp.N(Sq.subs(Pi, Pi_star), 30)
Sigma_star = sp.N(Sigma0.subs(Pi, Pi_star), 30)
That_star = sp.N(That.subs(Pi, Pi_star), 30)

print("g_minus^F1 =", sp.N(gminus, 30))
print("Pi_*       =", Pi_star)
print("g(Pi_*)    =", g_star)
print("R_q(Pi_*)  =", Rq_star)
print("S_q(Pi_*)  =", Sq_star)
print("Sigma0(Pi_*) =", Sigma_star)
print("That(Pi_*)   =", That_star)

if abs(float(g_star - sp.N(gminus,30))) > 1e-12:
    raise AssertionError("Pi_* does not solve the compensation equation accurately enough.")

banner("STAGE 125 LEDGER")
print("Self-consistent Family-1 mouth branch:")
print("  Pi = Sigma0 * [1 - R_q(Pi) S_q(Pi)]")
print("  Sigma0(Pi) = Pi / (1 - R_q(Pi) S_q(Pi))")
print("  That(Pi)   = sqrt(9 Sigma0(Pi) / 20)")
print()
print("Canonical point:")
print(f"  Pi_*       = {Pi_star}")
print(f"  Sigma0_*   = {Sigma_star}")
print(f"  That_*     = {That_star}")
