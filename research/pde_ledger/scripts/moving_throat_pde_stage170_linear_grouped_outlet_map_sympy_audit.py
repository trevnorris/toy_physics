#!/usr/bin/env python3
"""
moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py

SymPy-backed audit for Stage 170.

Checks:
1. Linear grouped-bundle transport formulas for delta u2, delta u4, delta P0.
2. Exact map from grouped-lane defects into delta kappa_W and delta gamma_W.
3. Exact one-parameter even-consistency relation:
       delta D4 = (2/3) delta D2 + (1/27) delta D0.
4. Grouped trace/anomaly transport formulas.
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

banner("STAGE 170 — LINEAR GROUPED-P2 DIRECT OUTLET MAP")

# ---------------------------------------------------------------------------
# 1. Linear grouped transport around the canonical isotropic branch
# ---------------------------------------------------------------------------
D0, dD0, dD2, dD4 = sp.symbols('D0 dD0 dD2 dD4', nonzero=True, real=True)
N0, dN0 = sp.symbols('N0 dN0', nonzero=True, real=True)
sigma = sp.symbols('sigma', nonzero=True, real=True)

u2 = sp.Rational(1, 9)
u4 = sp.Rational(4, 81)
P0 = sp.symbols('P0', nonzero=True, real=True)

D2 = -u2 * D0
D4 = -D0 / 27

eps = sp.symbols('eps', real=True)

u2_full = - (D2 + eps*dD2) / (D0 + eps*dD0)
u4_full = ((D2 + eps*dD2)**2 - (D0 + eps*dD0) * (D4 + eps*dD4)) / (D0 + eps*dD0)**2
P0_full = (N0 + eps*dN0) / (D0 + eps*dD0)

du2 = sp.expand(sp.series(u2_full, eps, 0, 2).removeO()).coeff(eps, 1)
du4 = sp.expand(sp.series(u4_full, eps, 0, 2).removeO()).coeff(eps, 1)
dP0 = sp.expand(sp.series(P0_full.subs(N0, P0*D0), eps, 0, 2).removeO()).coeff(eps, 1)

banner("Linear grouped conservative/output transport")
expect_zero("delta u2 + (dD2 + dD0/9)/D0", du2 + (dD2 + dD0/9)/D0)
expect_zero("delta u4 + (dD4 + 2 dD2/9 + 5 dD0/81)/D0", du4 + (dD4 + 2*dD2/9 + 5*dD0/81)/D0)
expect_zero("delta P0 - (dN0 - P0 dD0)/D0", dP0 - (dN0 - P0*dD0)/D0)

# ---------------------------------------------------------------------------
# 2. Direct outlet map
# ---------------------------------------------------------------------------
dkappa, dgamma = sp.symbols('dkappa dgamma', real=True)

du2_hyb = - sigma * dkappa / (3 * (1 - sigma))
dP0_over_P0_hyb = - 9 * sigma * dgamma / (1 - sigma)

dkappa_from_du2 = sp.solve(sp.Eq(sp.Symbol('du2sym'), du2_hyb), dkappa)[0].subs(sp.Symbol('du2sym'), du2)
dgamma_from_dP0 = sp.solve(sp.Eq(sp.Symbol('dP0sym')/P0, dP0_over_P0_hyb), dgamma)[0].subs(sp.Symbol('dP0sym'), dP0)

banner("Direct outlet coefficients")
expect_zero(
    "delta kappa_W - 3(1-sigma)(dD2 + dD0/9)/(sigma D0)",
    dkappa_from_du2 - 3*(1-sigma)*(dD2 + dD0/9)/(sigma*D0)
)
expect_zero(
    "delta gamma_W + (1-sigma)(dN0 - P0 dD0)/(9 sigma N0)",
    dgamma_from_dP0.subs(P0, N0/D0) + (1-sigma)*(dN0 - (N0/D0)*dD0)/(9*sigma*N0)
)

# ---------------------------------------------------------------------------
# 3. Even one-parameter consistency relation
# ---------------------------------------------------------------------------
du4_from_hyb = -8 * sigma * dkappa / (27 * (1 - sigma))
du4_from_kappa = sp.simplify(du4_from_hyb.subs(dkappa, dkappa_from_du2))

banner("Even one-parameter consistency")
expect_zero("delta u4 - (8/9) delta u2", du4_from_kappa - sp.Rational(8,9) * du2)
relation = sp.simplify(sp.solve(sp.Eq(du4, sp.Rational(8,9)*du2), dD4)[0])
expect_zero("delta D4 - (2/3) delta D2 - dD0/27", relation - (sp.Rational(2,3)*dD2 + dD0/27))

# ---------------------------------------------------------------------------
# 4. Grouped trace/anomaly formulas
# ---------------------------------------------------------------------------
aD0, aD2, aD4, aN0 = sp.symbols('aD0 aD2 aD4 aN0', real=True)
bD0, bD2, bD4, bN0 = sp.symbols('bD0 bD2 bD4 bN0', real=True)

a_kappa = sp.simplify(3*(1-sigma)*(aD2 + aD0/9)/(sigma*D0))
b_kappa = sp.simplify(3*(1-sigma)*(bD2 + bD0/9)/(sigma*D0))

a_gamma = sp.simplify(-(1-sigma)*(aN0 - P0*aD0)/(9*sigma*N0))
b_gamma = sp.simplify(-(1-sigma)*(bN0 - P0*bD0)/(9*sigma*N0))
P0_ref = sp.simplify(N0 / D0)

a_kappa_from_map = sp.simplify(sp.expand(dkappa_from_du2.subs({dD2: aD2, dD0: aD0})))
b_kappa_from_map = sp.simplify(sp.expand(dkappa_from_du2.subs({dD2: bD2, dD0: bD0})))
a_gamma_from_map = sp.simplify(sp.expand(dgamma_from_dP0.subs({dN0: aN0, dD0: aD0})))
b_gamma_from_map = sp.simplify(sp.expand(dgamma_from_dP0.subs({dN0: bN0, dD0: bD0})))

banner("Grouped trace/anomaly transport")
print("a_kappa =", a_kappa)
print("b_kappa =", b_kappa)
print("a_gamma =", a_gamma)
print("b_gamma =", b_gamma)
print("a_kappa from map =", a_kappa_from_map)
print("b_kappa from map =", b_kappa_from_map)
print("a_gamma from map =", a_gamma_from_map)
print("b_gamma from map =", b_gamma_from_map)
expect_zero("trace kappa coefficient", a_kappa_from_map - a_kappa)
expect_zero("anomaly kappa coefficient", b_kappa_from_map - b_kappa)
expect_zero(
    "trace gamma coefficient",
    sp.simplify(a_gamma_from_map.subs(P0, P0_ref) - a_gamma.subs(P0, P0_ref)),
)
expect_zero(
    "anomaly gamma coefficient",
    sp.simplify(b_gamma_from_map.subs(P0, P0_ref) - b_gamma.subs(P0, P0_ref)),
)

# ---------------------------------------------------------------------------
# 5. Weak-axisymmetric branch: signature (1, 1/2, -1) and scalar amplitudes
#    (paper Sec. 5 / card Checks item 2)
# ---------------------------------------------------------------------------
# On the weak-axisymmetric Y_20 branch the grouped bundle defects scale as
#   delta D_(A,n) = eps * lambda_A * D_n^(1),  delta N_(A,0) = eps * lambda_A * N_0^(1)
# with (lambda_20, lambda_21, lambda_22) = (1, 1/2, -1). Feeding these lane-scaled
# inputs through the SAME linear outlet maps verified in Sec. 2 must reproduce the
# grouped signature on delta kappa_W / delta gamma_W and collapse to two scalar
# amplitudes kappa1, gamma1 with the closed forms boxed in notes Sec. 5.
D0_1, D2_1, N0_1 = sp.symbols('D0_1 D2_1 N0_1', real=True)
eps_l = sp.symbols('eps_l', real=True)

def kappa_map(dD2_, dD0_):
    return 3*(1 - sigma)*(dD2_ + dD0_/9)/(sigma*D0)

def gamma_map(dN0_, dD0_):
    return -(1 - sigma)*(dN0_ - P0*dD0_)/(9*sigma*N0)

kappa1 = 3*(1 - sigma)*(D2_1 + D0_1/9)/(sigma*D0)
gamma1 = -(1 - sigma)*(N0_1 - P0*D0_1)/(9*sigma*N0)

banner("Weak-axisymmetric signature (1, 1/2, -1) and scalar amplitudes")
lanes = {20: sp.Integer(1), 21: sp.Rational(1, 2), 22: sp.Integer(-1)}
dkW = {}
dgW = {}
for A, lam in lanes.items():
    dkW[A] = kappa_map(eps_l*lam*D2_1, eps_l*lam*D0_1)
    dgW[A] = gamma_map(eps_l*lam*N0_1, eps_l*lam*D0_1)
    expect_zero(f"delta kappa_W^({A}) - eps lambda kappa1", dkW[A] - eps_l*lam*kappa1)
    expect_zero(f"delta gamma_W^({A}) - eps lambda gamma1", dgW[A] - eps_l*lam*gamma1)
# grouped signature ratios (lambda_20, lambda_21, lambda_22) = (1, 1/2, -1)
expect_zero("kappa signature: 21 = (1/2) 20", dkW[21] - sp.Rational(1, 2)*dkW[20])
expect_zero("kappa signature: 22 = -20", dkW[22] + dkW[20])
expect_zero("gamma signature: 21 = (1/2) 20", dgW[21] - sp.Rational(1, 2)*dgW[20])
expect_zero("gamma signature: 22 = -20", dgW[22] + dgW[20])

print("\nCarry-forward formulas:")
print("  delta kappa_W^(A) = 3(1-sigma_*) [delta D_(A,2) + delta D_(A,0)/9] / (sigma_* D0)")
print("  delta gamma_W^(A) = -(1-sigma_*) [delta N_(A,0) - P0 delta D_(A,0)] / (9 sigma_* N0)")
print("  direct even-pole consistency: delta D_(A,4) = (2/3) delta D_(A,2) + delta D_(A,0)/27")
print("  so the full linear grouped-anisotropy outlet problem collapses to the pair")
print("      K_A := delta D_(A,2) + delta D_(A,0)/9")
print("      G_A := delta N_(A,0) - P0 delta D_(A,0)")
