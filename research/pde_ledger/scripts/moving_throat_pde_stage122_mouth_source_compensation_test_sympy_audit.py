#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(s):
    print("\n" + "="*88)
    print(s)
    print("="*88)

def expect_zero(name, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

def expect_nonzero(name, expr):
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr == 0:
        raise AssertionError(f"{name} is unexpectedly zero")

banner("STAGE 122 — NATURAL MOUTH-SOURCE BRANCH VS COMPENSATION FAMILY")

R = sp.Rational(37,20)
rF = sp.sqrt(12*R**2/sp.pi**2 - 1)
gminus = sp.simplify(rF - sp.sqrt(1+rF**2)/2)
gplus  = sp.simplify(rF + sp.sqrt(1+rF**2)/2)

print("r_F1 =", rF)
print("g_minus =", gminus)
print("g_plus  =", gplus)

# natural equal-normalized source branch
g_nat = sp.Integer(1)

comp_def = sp.simplify(1 + rF**2 - 4*(g_nat - rF)**2)
print("compensation defect on natural branch =", comp_def)
print("numeric defect =", sp.N(comp_def, 20))

print("delta g_minus =", sp.N(g_nat - gminus, 20))
print("delta g_plus  =", sp.N(gplus - g_nat, 20))

# Independent traction ratio from stage 119 boxed law  g = C / T_m,
# with C = sqrt(2 Z_q K_s)/(J_s c_s sqrt(mu_0 L_W)) the SAME background-fixed
# constant on every branch. Then T_m^(±)/T_m^nat = (C/g_±)/(C/g_nat) = g_nat/g_±.
# C is carried symbolically so its cancellation is verified, not assumed; g_nat
# is the equal-normalized natural-branch ansatz value (line 34), not 1/g_±.
C = sp.symbols("C", positive=True)
Tm_nat   = C / g_nat
Tm_minus = C / gminus
Tm_plus  = C / gplus
T_ratio_minus = sp.simplify(Tm_minus / Tm_nat)
T_ratio_plus  = sp.simplify(Tm_plus  / Tm_nat)
print("T_m(-)/T_m(nat) =", T_ratio_minus)
print("T_m(+)/T_m(nat) =", T_ratio_plus)
print("numeric T ratio (-) =", sp.N(T_ratio_minus, 20))
print("numeric T ratio (+) =", sp.N(T_ratio_plus, 20))

# exact closed forms
gminus_exact = sp.simplify((2*sp.sqrt(4107 - 100*sp.pi**2) - 37*sp.sqrt(3)) / (20*sp.pi))
gplus_exact = sp.simplify((2*sp.sqrt(4107 - 100*sp.pi**2) + 37*sp.sqrt(3)) / (20*sp.pi))
expect_zero("gminus exact form", gminus - gminus_exact)
expect_zero("gplus exact form", gplus - gplus_exact)

expect_zero("compensation quadratic at gminus", 1 + rF**2 - 4*(gminus - rF)**2)
expect_zero("compensation quadratic at gplus",  1 + rF**2 - 4*(gplus  - rF)**2)

defect_exact = (-12321 + 80*sp.pi*sp.sqrt(4107 - 100*sp.pi**2)) / (100*sp.pi**2)
expect_zero("defect closed form", comp_def - defect_exact)

expect_nonzero("natural off compensation", comp_def)

expect_zero("traction ratio (-) identity", T_ratio_minus - 1/gminus)
expect_zero("traction ratio (+) identity", T_ratio_plus  - 1/gplus)
