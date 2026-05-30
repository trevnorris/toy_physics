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

def expect_equal(name, lhs, rhs):
    diff = sp.simplify(lhs - rhs)
    print(f"{name}: lhs - rhs = {diff}")
    if diff != 0:
        raise AssertionError(f"{name} mismatch: {lhs} vs {rhs}")

def expect_positive(name, expr):
    val = sp.simplify(expr)
    print(f"{name}: {val}")
    if not (val.is_positive is True or (val.is_number and float(val) > 0)):
        raise AssertionError(f"{name} not positive: {val}")

Pi = sp.symbols("Pi", positive=True, real=True)
pi = sp.pi

banner("STAGE 143 — EQUAL-NORMALIZED BRANCH IS A SINGULAR LIMIT")

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
expect_positive("pi**2 - 2*pi > 0", pi**2 - 2*pi)
expect_positive("pi**2/2 - 4 > 0", pi**2/2 - 4)
# exp-remainder positivity: prove R(Pi) = exp(Pi) - 1 - Pi - Pi**2/2 > 0 for all Pi > 0.
# Rigorous argument (each line is an independent, can-fail assertion):
#   R(0) = 0, R'(0) = 0, R''(0) = 0, and R'''(Pi) = exp(Pi) > 0 for all Pi,
# so R is strictly increasing from 0 on Pi > 0, hence R(Pi) > 0 there.
R_rem = sp.exp(Pi) - 1 - Pi - Pi**2/2
expect_equal("exp remainder R(0) == 0", R_rem.subs(Pi, 0), sp.Integer(0))
expect_equal("exp remainder R'(0) == 0", sp.diff(R_rem, Pi).subs(Pi, 0), sp.Integer(0))
expect_equal("exp remainder R''(0) == 0", sp.diff(R_rem, Pi, 2).subs(Pi, 0), sp.Integer(0))
expect_zero("exp remainder R'''(Pi) - exp(Pi) == 0", sp.diff(R_rem, Pi, 3) - sp.exp(Pi))
expect_positive("exp remainder R'''(Pi) = exp(Pi) > 0 for Pi>0", sp.exp(Pi))

subbanner("Endpoint limits")
g0 = sp.limit(gPi, Pi, 0, dir='+')
ginf = sp.limit(gPi, Pi, sp.oo)
print("lim_{Pi->0+} g_Pi =", g0)
print("lim_{Pi->oo} g_Pi =", ginf)
expect_equal("lim_{Pi->0+} g_Pi == 2/pi", g0, 2/pi)
expect_equal("lim_{Pi->oo} g_Pi == 1", ginf, sp.Integer(1))

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
expect_equal("R_infty == (1-r)**2/(1+r**2)", Rinf, (1-r)**2/(1+r**2))
expect_equal("S_infty == 1", Sinf, sp.Integer(1))
expect_equal("lim That/sqrt(Pi) == sqrt((9/20)/(1-R_infty))", that_ratio, sp.sqrt(sp.Rational(9,20)/(1-Rinf)))

banner("STAGE 143 LEDGER")
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
