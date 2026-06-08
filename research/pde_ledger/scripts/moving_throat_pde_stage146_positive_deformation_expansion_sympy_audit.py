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

# Symbolic verification that the closed-form moments equal the canonical
# integrals against the kernel functions. This anchors the formulas before
# the numeric Pi_* root-finder runs.
gPi_direct = sp.integrate(Sigma*sp.cos(sp.pi*x/2), (x, 0, 1))
if gPi_direct.has(sp.Integral):
    for pp in [sp.Rational(7,10), sp.Rational(11,10), sp.Rational(17,10), sp.Rational(23,10)]:
        diff_g = sp.N(sp.integrate((Sigma*sp.cos(sp.pi*x/2)).subs(Pi, pp), (x,0,1)) - gPi.subs(Pi, pp), 30)
        print(f"g(Pi) numeric sample Pi={pp}: diff={diff_g}")
        if abs(float(diff_g)) > 1e-25:
            raise AssertionError(f"g(Pi) numeric mismatch at Pi={pp}: {diff_g}")
    print("PASS: g(Pi) direct-formula (numeric fallback at 4 samples)")
else:
    expect_zero("g(Pi) direct-formula", sp.simplify(gPi_direct - gPi))

Sq_direct = sp.integrate(Sigma*Kq, (x, 0, 1))
if Sq_direct.has(sp.Integral):
    for pp in [sp.Rational(7,10), sp.Rational(11,10), sp.Rational(17,10), sp.Rational(23,10)]:
        diff_s = sp.N(sp.integrate((Sigma*Kq).subs(Pi, pp), (x,0,1)) - Sformula.subs(Pi, pp), 30)
        print(f"S_q(Pi) numeric sample Pi={pp}: diff={diff_s}")
        if abs(float(diff_s)) > 1e-25:
            raise AssertionError(f"S_q(Pi) numeric mismatch at Pi={pp}: {diff_s}")
    print("PASS: S_q(Pi) direct-formula (numeric fallback at 4 samples)")
else:
    expect_zero("S_q(Pi) direct-formula", sp.simplify(Sq_direct - Sformula))

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
Pi_star = sp.nsolve(gPi - gminus, 1.5, prec=50)
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

# Convex-family affine moment laws tested by direct quadrature of the
# assembled profile against closed-form intercepts g_* and S_*.
# Use a concrete positive normalized test profile varsigma_test(x) =
# 6 x (1 - x), positive on (0,1), with integral 1 on [0,1].
varsigma_test = 6*x*(1-x)
Sigma_eps = (1-eps)*Sigma.subs(Pi, Pi_star) + eps*varsigma_test
def quad_moment(expr):
    return sp.N(sp.Integral(expr, (x, 0, 1)).evalf(50), 50)

g_star_closed = sp.N(gPi.subs(Pi, Pi_star), 50)
S_star_closed = sp.N(Sformula.subs(Pi, Pi_star), 50)
gbar_v = quad_moment(varsigma_test*sp.cos(sp.pi*x/2))
Sbar_v = quad_moment(varsigma_test*Kq)
g_slope = sp.N(gbar_v - g_star_closed, 50)
S_slope = sp.N(Sbar_v - S_star_closed, 50)
TOL = sp.Float("1e-25", 50)
SLOPE_TOL = sp.Float("1e-3", 50)
print(f"g_eps convex affine moment law nonzero slope |gbar_v - g_*| = {sp.Abs(g_slope)}")
print(f"S_eps convex affine moment law nonzero slope |Sbar_v - S_*| = {sp.Abs(S_slope)}")
if sp.Abs(g_slope) <= SLOPE_TOL:
    raise AssertionError(f"g_eps affine slope is vacuous: |gbar_v - g_*|={sp.Abs(g_slope)} <= 1e-3")
if sp.Abs(S_slope) <= SLOPE_TOL:
    raise AssertionError(f"S_eps affine slope is vacuous: |Sbar_v - S_*|={sp.Abs(S_slope)} <= 1e-3")
for eps_val, label in [(sp.Rational(1, 10), "eps=1/10"), (sp.Rational(1, 2), "eps=1/2")]:
    Sigma_eps_sample = Sigma_eps.subs(eps, eps_val)
    gbar_eps = quad_moment(Sigma_eps_sample*sp.cos(sp.pi*x/2))
    Sbar_eps = quad_moment(Sigma_eps_sample*Kq)
    g_res = sp.N(gbar_eps - (g_star_closed + eps_val*(gbar_v - g_star_closed)), 50)
    S_res = sp.N(Sbar_eps - (S_star_closed + eps_val*(Sbar_v - S_star_closed)), 50)
    print(f"g_eps convex affine moment law via direct quadrature at {label}: residual = {g_res}")
    print(f"S_eps convex affine moment law via direct quadrature at {label}: residual = {S_res}")
    if sp.Abs(g_res) >= TOL:
        raise AssertionError(f"g_eps convex affine moment law via direct quadrature fails at {label}: residual={g_res} >= 1e-25")
    if sp.Abs(S_res) >= TOL:
        raise AssertionError(f"S_eps convex affine moment law via direct quadrature fails at {label}: residual={S_res} >= 1e-25")
print("PASS: g_eps convex affine moment law via direct quadrature with closed-form g_* intercept and nonzero slope")
print("PASS: S_eps convex affine moment law via direct quadrature with closed-form S_* intercept and nonzero slope")

print("\nStage 146 complete.")
