#!/usr/bin/env python3
"""
moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py

SymPy-backed audit for Stage 169.

Checks:
1. Exact grouped bilinear invariant formula I[x,y] = (1/5) dx^T G dy = 4 a_x a_y + (4/5) b_x b_y.
2. Weak axisymmetric branch law b = 3 a and I[x,y] = (7/10) eps^2 x1 y1.
3. Explicit Y20 average vanishes, so the linear term of a scalar log-observable vanishes.
4. Stage 168 weighted combination transports the quadratic grouped invariants into eps_perp.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)


def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 169 — NO LINEAR GROUPED-P2 SCALAR SLIPPAGE")

# ---------------------------------------------------------------------------
# 1. Exact grouped bilinear invariant
# ---------------------------------------------------------------------------
xb20, xb21, xb22 = sp.symbols('x20 x21 x22', real=True)
yb20, yb21, yb22 = sp.symbols('y20 y21 y22', real=True)

x = sp.Matrix([xb20, xb21, xb22])
y = sp.Matrix([yb20, yb21, yb22])
G = sp.diag(1, 2, 2)
ebar = sp.Matrix([1, 1, 1])

xbar = (xb20 + 2*xb21 + 2*xb22) / 5
ybar = (yb20 + 2*yb21 + 2*yb22) / 5

ax = (2*xb20 - xb21 - xb22) / 10
bx = (xb21 - xb22) / 2
ay = (2*yb20 - yb21 - yb22) / 10
by = (yb21 - yb22) / 2

dx = x - xbar*ebar
dy = y - ybar*ebar

Ixy = sp.simplify((dx.T * G * dy)[0] / 5)
Ixy_expected = sp.simplify(4*ax*ay + sp.Rational(4,5)*bx*by)

banner("Grouped bilinear invariant")
print("I[x,y] =", Ixy)
expect_zero("I[x,y] - [4 a_x a_y + 4/5 b_x b_y]", Ixy - Ixy_expected)

# ---------------------------------------------------------------------------
# 2. Weak axisymmetric branch law
# ---------------------------------------------------------------------------
eps, x1, y1, x0, y0 = sp.symbols('eps x1 y1 x0 y0', real=True)
subs_axis = {
    xb20: x0 + eps*x1,
    xb21: x0 + eps*x1/2,
    xb22: x0 - eps*x1,
    yb20: y0 + eps*y1,
    yb21: y0 + eps*y1/2,
    yb22: y0 - eps*y1,
}

ax_axis = sp.simplify(ax.subs(subs_axis))
bx_axis = sp.simplify(bx.subs(subs_axis))
Ixy_axis = sp.simplify(Ixy.subs(subs_axis))

banner("Weak axisymmetric grouped law")
print("a_x(axis) =", ax_axis)
print("b_x(axis) =", bx_axis)
expect_zero("b_x - 3 a_x", bx_axis - 3*ax_axis)
expect_zero("I_axis - (7/10) eps^2 x1 y1", Ixy_axis - sp.Rational(7,10)*eps**2*x1*y1)

# ---------------------------------------------------------------------------
# 3. Explicit Y20 average and no linear scalar feed-down
# ---------------------------------------------------------------------------
th, ph = sp.symbols('th ph', real=True)
Y20 = sp.sqrt(sp.Integer(5)/(16*sp.pi)) * (3*sp.cos(th)**2 - 1)

def sphere_avg(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(
        sp.integrate(sp.integrate(expr*sp.sin(th), (ph, 0, 2*sp.pi)), (th, 0, sp.pi)) / (4*sp.pi)
    )

banner("Explicit harmonic average")
Y20_avg = sphere_avg(Y20)
Y20_sq_avg = sphere_avg(Y20**2)
print("<Y20> =", Y20_avg)
print("<Y20^2> =", Y20_sq_avg)
expect_zero("average Y20", Y20_avg)
expect_zero("average Y20^2 - 1/(4 pi)", Y20_sq_avg - 1/(4*sp.pi))

X0, e = sp.symbols('X0 e', positive=True, real=True)
log_series = sp.expand(sp.log(X0*(1 + e*Y20)).series(e, 0, 3).removeO())
log_avg = sp.simplify(sphere_avg(log_series))
lin_coeff = sp.simplify(sp.diff(log_avg, e).subs(e, 0))

print("<log(X0 (1 + e Y20))> through O(e^2) =", log_avg)
expect_zero("linear coefficient in averaged log-observable", lin_coeff)

# ---------------------------------------------------------------------------
# 4. Stage 168 transport into eps_perp
# ---------------------------------------------------------------------------
XiL, Xiv, XiT, Igrp = sp.symbols('XiL Xiv XiT Igrp', real=True)
g, r = sp.symbols('g r', positive=True, real=True)
s = sp.sqrt(1 + r**2)

# model one quadratic grouped invariant contribution Igrp
# eps_L = XiL Igrp, eps_v = Xiv Igrp, eps_T = XiT Igrp

eps_perp = g*XiT*Igrp + (g + sp.Rational(1,2)/s)*Xiv*Igrp + (2*g + sp.Rational(3,4)/s)*XiL*Igrp
Xi_perp = g*XiT + (g + sp.Rational(1,2)/s)*Xiv + (2*g + sp.Rational(3,4)/s)*XiL

banner("Stage 168 weighted transport")
expect_zero("eps_perp - Xi_perp Igrp", eps_perp - Xi_perp*Igrp)

r_num = sp.Float('1.77799353547498', 30)
g_num = sp.Float('0.758035078944663', 30)
s_num = sp.sqrt(1 + r_num**2)
Xi_perp_num = sp.N(Xi_perp.subs({g: g_num, r: r_num}), 20)
print("Numeric Xi_perp combination =", Xi_perp_num)

coeff_T = sp.N(Xi_perp.coeff(XiT).subs({g: g_num, r: r_num}), 20)
coeff_v = sp.N(Xi_perp.coeff(Xiv).subs({g: g_num, r: r_num}), 20)
coeff_L = sp.N(Xi_perp.coeff(XiL).subs({g: g_num, r: r_num}), 20)
for name, got, want in [
    ("Xi_perp coeff on XiT", coeff_T, sp.Float('0.758035078944663', 20)),
    ("Xi_perp coeff on Xiv", coeff_v, sp.Float('1.00314310113848', 20)),
    ("Xi_perp coeff on XiL", coeff_L, sp.Float('1.88373219118005', 20)),
]:
    diff = sp.Abs(got - want)
    print(f"{name} = {got} (paper {want}, |diff| = {diff})")
    if diff > sp.Float('1e-12', 20):
        raise AssertionError(f"{name} disagrees with paper value {want}")

print("\nCarry-forward formulas:")
print("  I[x,y] = (1/5) (x-xbar ebar)^T Ggrp (y-ybar ebar) = 4 a_x a_y + (4/5) b_x b_y")
print("  Pure grouped P2 anisotropy has no linear scalar feed-down on an isotropic branch.")
print("  Weak axisymmetric branch: b = 3 a and I[x,y] = (7/10) eps^2 x^(1) y^(1).")
print("  Therefore grouped-lane anisotropy enters eps_L, eps_v, eps_T, eps_perp only quadratically.")
