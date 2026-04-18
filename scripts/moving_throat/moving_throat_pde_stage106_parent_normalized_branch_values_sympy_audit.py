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

banner("STAGE 106 — PARENT-NORMALIZED BRANCH VALUES")

# Symbols
a, L, ell, mu0, Zq, m, rho, q, hbar, cs = sp.symbols(
    'a L ell mu0 Z_q m rho q hbar c_s', positive=True, real=True
)
r, g = sp.symbols('r g', positive=True, real=True)

# Healing-locked shell formulas and D/N tube with L_W = L
Ks = 3*sp.pi*a**2*hbar**2/(5*m*rho*ell)
Kq = Zq/mu0 * sp.pi**2 * cs**2/(4*L**2)
Js = 4*sp.pi*a**2*ell/3
lam = -(8*sp.sqrt(2)/3) * q * sp.symbols('v_w0', real=True) * a**2 * ell * sp.sqrt(L)  # only for reference

# Parent-normalized Xi_v:
vw = sp.symbols('v_w0', real=True)
r_from_parent = sp.simplify(lam.subs({sp.symbols('v_w0'): vw})/sp.sqrt(Ks*Kq))
Xi_v = sp.symbols('Xi_v', real=True)
Xi_v_def = sp.simplify(q*sp.sqrt(mu0*m*rho)*a*L**sp.Rational(3,2)*ell**sp.Rational(3,2)*vw/(hbar*sp.sqrt(Zq)*cs))
Xi_v_expr = sp.simplify(Xi_v_def.subs(vw, sp.solve(sp.Eq(r, r_from_parent), vw)[0]))
print("Xi_v(r) =", Xi_v_expr)

# Parent-normalized Xi_T:
Tm = sp.symbols('T_m', positive=True, real=True)
g_from_parent = sp.simplify(sp.sqrt(2*Zq*Ks)/(Tm*Js*cs*sp.sqrt(mu0*L))).subs({cs: hbar/(2*m*ell)})
Xi_T = sp.symbols('Xi_T', positive=True, real=True)
Xi_T_def = sp.simplify(Tm*a*sp.sqrt(mu0*rho*L*ell)/(sp.sqrt(Zq*m)))
Xi_T_expr = sp.simplify(Xi_T_def.subs(Tm, sp.solve(sp.Eq(g, g_from_parent), Tm)[0]))
print("Xi_T(g) =", Xi_T_expr)

expect_zero("Xi_v law", Xi_v_expr + 3*sp.sqrt(30)*sp.pi**sp.Rational(3,2)*r/160)
expect_zero("Xi_T law", Xi_T_expr - 3*sp.sqrt(30)/(10*sp.sqrt(sp.pi)*g))

# Family-1 numerical values
R = sp.Rational(37,20)
rF = sp.sqrt(12*R**2/sp.pi**2 - 1)
gminus = sp.simplify(rF - sp.sqrt(1+rF**2)/2)
gplus  = sp.simplify(rF + sp.sqrt(1+rF**2)/2)
Xi_v_F1 = sp.simplify((-3*sp.sqrt(30)*sp.pi**sp.Rational(3,2)/160) * rF)
Xi_T_nat = sp.simplify(3*sp.sqrt(30)/(10*sp.sqrt(sp.pi)))
Xi_T_minus = sp.simplify(Xi_T_nat/gminus)
Xi_T_plus  = sp.simplify(Xi_T_nat/gplus)

print("Xi_v(F1) =", Xi_v_F1, "≈", sp.N(Xi_v_F1, 20))
print("Xi_T(nat) =", Xi_T_nat, "≈", sp.N(Xi_T_nat, 20))
print("Xi_T(-) =", Xi_T_minus, "≈", sp.N(Xi_T_minus, 20))
print("Xi_T(+) =", Xi_T_plus, "≈", sp.N(Xi_T_plus, 20))
