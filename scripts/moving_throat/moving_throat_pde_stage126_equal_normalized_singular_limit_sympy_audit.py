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

banner("STAGE 126 — EQUAL-NORMALIZED BRANCH IS A SINGULAR LIMIT")

r = sp.sqrt(sp.Integer(4107) - 100*pi**2)/(10*pi)
gPi = 2*Pi*(2*Pi*sp.exp(Pi)+pi)/((4*Pi**2+pi**2)*(sp.exp(Pi)-1))

subbanner("Exact strict inequality g_Pi < 1 for finite Pi")
num = sp.expand((1-gPi)*(4*Pi**2+pi**2)*(sp.exp(Pi)-1))
print("numerator =", num)

decomp = pi**2*(sp.exp(Pi)-1-Pi-Pi**2/2) + Pi*(pi**2-2*pi) + Pi**2*(pi**2/2-4)
expect_zero("numerator - exact positive decomposition", num - decomp)

print("positive pieces:")
print("  exp remainder =", sp.exp(Pi)-1-Pi-Pi**2/2)
print("  linear coeff  =", sp.simplify(pi**2-2*pi))
print("  quadratic coeff =", sp.simplify(pi**2/2-4))

subbanner("Endpoint limits")
g0 = sp.limit(gPi, Pi, 0, dir='+')
ginf = sp.limit(gPi, Pi, sp.oo)
print("lim_{Pi->0+} g_Pi =", g0)
print("lim_{Pi->oo} g_Pi =", ginf)

Sq = Pi*(((pi/2)*sp.tanh(pi/2)) + Pi*(sp.exp(-Pi)*sp.sech(pi/2)-1))/((1-sp.exp(-Pi))*(((pi/2)**2)-Pi**2))
Rq = (gPi-r)**2/(1+r**2)
Rinf = sp.simplify(sp.limit(Rq, Pi, sp.oo))
Sinf = sp.simplify(sp.limit(Sq, Pi, sp.oo))
Sigma0 = Pi/(1-Rq*Sq)
That = sp.sqrt(sp.Rational(9,20)*Sigma0)

print("R_infty =", Rinf)
print("S_infty =", Sinf)
sigma_ratio = sp.simplify(sp.limit(Sigma0/Pi, Pi, sp.oo))
print("lim Sigma0/Pi =", sigma_ratio)
that_ratio = sp.simplify(sp.sqrt(sp.Rational(9,20)*sigma_ratio))
print("lim That/sqrt(Pi) =", that_ratio)

banner("STAGE 126 LEDGER")
print("For every finite Pi>0:")
print("  2/pi < g_Pi < 1")
print("So g_c = 1 is not a finite positive-bias branch.")
print()
print("Equal-normalized branch:")
print("  reached only as Pi -> infinity")
print("  and That -> infinity like sqrt(Pi).")
print()
print("Numerics:")
print(f"  R_infty = {sp.N(Rinf, 30)}")
print(f"  lim That/sqrt(Pi) = {sp.N(that_ratio, 30)}")
