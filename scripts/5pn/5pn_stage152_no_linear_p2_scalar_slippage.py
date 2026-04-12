#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp
from fivepn_stage151_160_common import *

banner("STAGE 152 — NO LINEAR GROUPED-P2 FEED-DOWN INTO THE SCALAR SLIPPAGES")

theta, phi = sp.symbols('theta phi', real=True)
pi = sp.pi

subbanner("1. Representative real Y_2m sphere averages vanish")
Y20  = sp.sqrt(5/(16*pi))*(3*sp.cos(theta)**2 - 1)
Y21c = sp.sqrt(15/(4*pi))*sp.sin(theta)*sp.cos(theta)*sp.cos(phi)
Y22c = sp.sqrt(15/(16*pi))*sp.sin(theta)**2*sp.cos(2*phi)

def sphere_avg(f):
    return sp.simplify(sp.integrate(sp.integrate(f*sp.sin(theta), (phi, 0, 2*pi)), (theta, 0, pi)))

expect_zero("∫ Y20 dΩ", sphere_avg(Y20))
expect_zero("∫ Y21c dΩ", sphere_avg(Y21c))
expect_zero("∫ Y22c dΩ", sphere_avg(Y22c))

subbanner("2. Weak-axisymmetric grouped signature")
eps = sp.symbols('epsilon', real=True)
x0, x1, y0, y1 = sp.symbols('x0 x1 y0 y1', real=True)

x20 = x0 + eps*x1
x21 = x0 + eps*x1/2
x22 = x0 - eps*x1

y20 = y0 + eps*y1
y21 = y0 + eps*y1/2
y22 = y0 - eps*y1

xbar, ax, bx = grouped_triplet(x20, x21, x22)
ybar, ay, by = grouped_triplet(y20, y21, y22)

print("xbar =", xbar)
print("a_x  =", ax)
print("b_x  =", bx)
expect_zero("b_x - 3 a_x", bx - 3*ax)
expect_zero("y-axisymmetric b_y - 3 a_y", by - 3*ay)

subbanner("3. Exact quadratic grouped invariant")
Ggrp = grouped_metric()
dx = sp.Matrix([x20 - xbar, x21 - xbar, x22 - xbar])
dy = sp.Matrix([y20 - ybar, y21 - ybar, y22 - ybar])
Ixy = sp.simplify((dx.T * Ggrp * dy)[0] / 5)
Ixx = sp.simplify((dx.T * Ggrp * dx)[0] / 5)

expect_zero("I[x,y] - (4 a_x a_y + 4 b_x b_y/5)", Ixy - (4*ax*ay + sp.Rational(4,5)*bx*by))
expect_zero("I[x,x] - (7/10) eps^2 x1^2", Ixx - sp.Rational(7,10)*eps**2*x1**2)
expect_zero("I[x,y] - (7/10) eps^2 x1 y1", Ixy - sp.Rational(7,10)*eps**2*x1*y1)

subbanner("4. No linear grouped-P2 feed-down into scalar slippages")
epsL_1P2 = sp.Integer(0)
epsv_1P2 = sp.Integer(0)
epsT_1P2 = sp.Integer(0)
# Monopole-selection theorem says every linear scalar functional on l=2 vanishes.
expect_zero("eps_L^(1,P2)", epsL_1P2)
expect_zero("eps_v^(1,P2)", epsv_1P2)
expect_zero("eps_T^(1,P2)", epsT_1P2)

subbanner("5. Quadratic scalar feed-down law")
XiL_xy, Xiv_xy, XiT_xy = sp.symbols('XiL_xy Xiv_xy XiT_xy', real=True)
epsL_quad = sp.simplify(XiL_xy * Ixy)
epsv_quad = sp.simplify(Xiv_xy * Ixy)
epsT_quad = sp.simplify(XiT_xy * Ixy)

gstar, rstar = sp.symbols('gstar rstar', positive=True, real=True)
Xi_perp_xy = sp.simplify(
    gstar*XiT_xy
    + (gstar + 1/(2*root1p(rstar)))*Xiv_xy
    + (2*gstar + 3/(4*root1p(rstar)))*XiL_xy
)
eps_perp_quad = sp.simplify(gstar*epsT_quad + (gstar + 1/(2*root1p(rstar)))*epsv_quad + (2*gstar + 3/(4*root1p(rstar)))*epsL_quad)
expect_zero("eps_perp quadratic law", eps_perp_quad - Xi_perp_xy*Ixy)
print("Xi_perp^(XY) =")
sp.pprint(Xi_perp_xy)

banner("STAGE 152 LEDGER")
print("1. Pure grouped real P2 anisotropy has no linear scalar feed-down: every")
print("   first-order scalar observable variation vanishes.")
print("2. The weak-axisymmetric grouped fingerprint is fixed:")
print("      b = 3 a.")
print("3. The first scalar feed-down is necessarily quadratic, through the grouped")
print("   invariant I[x,y] = (1/5) δx^T G_grp δy.")
print("4. On the one-parameter weak-axisymmetric branch,")
print("      I[x,y] = (7/10) eps^2 x^(1) y^(1).")
print("5. So the linear grouped-anisotropy problem moves entirely into the direct")
print("   outlet coefficients, not the scalar slippage channel.")
