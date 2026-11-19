"""Test whether a known axial pull could select the guide flow and its head."""

import sympy as sp

a, rho, ratio, D, C, sigma, ell, z = sp.symbols("a rho m D C sigma ell z", positive=True)
area = 4*sp.pi*a**3/3
pull = C+4*sp.pi*sigma*a**2
speed = D/(rho*area)
thrust_difference = rho*area*speed**2*(1-1/ratio)**2/2
rate_sq = sp.solve(thrust_difference-pull,D**2)[0]
head = sp.factor((rho*speed**2/(2*ratio**2)).subs(D**2,rate_sq))
print("mass-rate squared required by axial balance =",sp.factor(rate_sq))
print("required inlet total head above outlet ambient, pressure units =",head)
print("emission from the balanced normal load =",sp.factor((ell*z*thrust_difference).subs(D**2,rate_sq)))
head_ratio = sp.factor(head.subs(a,2*a)/head)
print("head at twice the radius / head at original radius =",head_ratio)
print("curvature-only and tension-only head ratios =",head_ratio.subs(sigma,0),head_ratio.subs(C,0))
assert sp.factor(head_ratio-(1+3*(4*sp.pi*sigma*a**2/pull))/8)==0
assert sp.factor(head*area*(ratio-1)**2-pull)==0

# For a fixed negative target B/L, solve the radius-dependent geometry that
# would be needed. This is a compatibility relation for independent data.
K = sp.symbols("K", positive=True)
large_root=(K+sp.sqrt(2*K-1))/(K-1)
gain_factor=(ratio**2+1)/(ratio-1)**2
print("widening area ratio needed for local thrust factor K>1 =",large_root)
assert sp.simplify(gain_factor.subs(ratio,large_root)-K)==0
